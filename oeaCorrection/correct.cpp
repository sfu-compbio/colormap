/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#include "correct.h"
#include "reader.h"
#include <pthread.h>

pthread_t				*_cor_threads = NULL;
pthread_mutex_t 		_cor_printLock;
ostringstream 			*_cor_outBuffer;
const double 			_cor_max_weight = std::numeric_limits<double>::infinity();
int 					MAX_MISMATCH = 1;
map<string, string> 	_cor_allSRs;

void* 					correctLrST(void *idp);
void 					correctSingleLr(int tid, read_t *rp, deque<sam_t> &i2pList, map<string, sam_t>& allMappings);
void 					getComponents(string &lr, deque<component_t> &lst);
void 					storeReadNIndex(deque<component_t> &compList, int tid);
void 					extractOEAs(deque<sam_t> &i2pList, map<string, sam_t> &allMappings, set<string> &readList, int bPos, int ePos, bool isForward);
void 					construtContigsFromOEA(set<string> &readList, deque<sam_t> &mapList, string fileName, map<string, string> &srList, int pivot, int tid);
string 					getQuerySubstr(sam_t &a, int rs, int re);
void 					updateCompPos(deque<component_t> &compList, int index, int offset);
string 					getMateId(string id);
void 					getAllSamEntries(string path, deque<sam_t> &lst, int pivot);

void initCorrect()
{
	// get all short reads
	getAllReads(SR_FILE, _cor_allSRs);

	_cor_threads = (pthread_t*) malloc(THREAD_COUNT * sizeof(pthread_t));
	_cor_outBuffer = new ostringstream[THREAD_COUNT];
}

void finalizeCorrect()
{
	free(_cor_threads);
	delete[] _cor_outBuffer;
}

void correctLrMT()
{
	int i;

	for (i = 0; i < THREAD_COUNT; i++)
	{
		_cor_outBuffer[i].str("");
		_cor_outBuffer[i].clear();
	}

	for (i = 0; i < THREAD_COUNT; i++)
		pthread_create(_cor_threads + i, NULL, correctLrST, THREAD_ID + i);

	for (i = 0; i < THREAD_COUNT; i++)
		pthread_join(_cor_threads[i], NULL);

	for (i = 0; i < THREAD_COUNT; i++)
		cout<< _cor_outBuffer[i].str();
}

void* correctLrST(void *idp)
{
	int tid;
	tid = *(int*)idp;

	read_t *rp;
	deque<sam_t> samList;
	map<string, sam_t> allMappings;
	while((rp=getNextReadsInfo(samList, allMappings)) != NULL)
	{
		correctSingleLr(tid, rp, samList, allMappings);
	}
}

void correctSingleLr(int tid, read_t *rp, deque<sam_t> &i2pList, map<string, sam_t>& allMappings)
{
	// deque<m5_t> nodeList;
	deque<component_t> compList;
	string currentLrId = rp->id;
	string &currentLrSeq = rp->seq;
	set<string> readList;
	deque<sam_t> mapList_left;
	deque<sam_t> mapList_right;
	string corrLR;
	string tmpStr;
	// int cntUncorrComp;
	// int checkCnt=0;
	// int counter=0;

	if(i2pList.size() > 0)
	{
		int i;
		getComponents(currentLrSeq, compList);
		
		int offset;
		int tmpPos;
		for(i=0; i<compList.size(); i++)
		{
			// if in the middle of two corrected regions
			if(compList[i].isCorrected == false && 
					i-1>=0 && compList[i-1].isCorrected == true &&
					i+1<compList.size() && compList[i+1].isCorrected == true)
			{
				storeReadNIndex(compList, tid);
				extractOEAs(i2pList, allMappings, readList, compList[i-1].rStart, compList[i-1].rEnd, true);
				construtContigsFromOEA(readList, mapList_left, "oea_left.fasta", _cor_allSRs, compList[i].rStart, tid);
				extractOEAs(i2pList, allMappings, readList, compList[i+1].rStart, compList[i+1].rEnd, false);
				construtContigsFromOEA(readList, mapList_right, "oea_right.fasta", _cor_allSRs, compList[i].rEnd, tid);
				// exit(0);

				// extend 
				if(mapList_left.size()>0 && mapList_right.size()>0) // both left and right
				{
					// if(mapList_left.size()>1)
					// 	checkCnt++;
					// if(mapList_right.size()>1)
					// 	checkCnt++;
					if(mapList_left[0].rEnd >= mapList_right[0].rStart) // overlap, that's great!!!
					{
						// cerr<< "case1" << endl;
						int posToExtend;
						posToExtend = min(mapList_left[0].rEnd, compList[i].rEnd);
						tmpStr = getQuerySubstr(mapList_left[0], compList[i].rStart, posToExtend);
						tmpStr += getQuerySubstr(mapList_right[0], posToExtend+1, compList[i].rEnd);
						// tmpStr is the corrected string of the un-corrected component
						compList[i-1].seq += tmpStr + compList[i+1].seq;
						// adjust positions
						tmpPos = compList[i-1].rStart + compList[i-1].seq.size() - 1;
						offset = tmpPos - compList[i+1].rEnd;
						compList[i-1].rEnd = tmpPos;
						updateCompPos(compList, i+2, offset);
						compList.erase(compList.begin()+i, compList.begin()+i+2);
						i--;
						/////////////////
						// cerr<< "======== new compList ========" << endl;
						// cerr<< offset << endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].seq;
						// }
						// cerr<< endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].rStart << " " << compList[z].rEnd << endl;
						// }
						/////////////////
					}
					else // no overlap, process left and right separately
					{
						// cerr<< "case2" << endl;
						// left
						tmpStr = getQuerySubstr(mapList_left[0], compList[i].rStart, mapList_left[0].rEnd);
						compList[i-1].seq += tmpStr;
						compList[i].seq = compList[i].seq.substr(mapList_left[0].rEnd-compList[i].rStart+1);
						// right
						tmpStr = getQuerySubstr(mapList_right[0], mapList_right[0].rStart, compList[i].rEnd);
						compList[i+1].seq = tmpStr + compList[i+1].seq;
						compList[i].seq = compList[i].seq.substr(0, compList[i].seq.size()-(compList[i].rEnd-mapList_right[0].rStart+1));
						// adjust positions
						compList[i-1].rEnd = compList[i-1].rStart + compList[i-1].seq.size() - 1;
						compList[i].rStart = compList[i-1].rEnd+1;
						compList[i].rEnd = compList[i].rStart + compList[i].seq.size() - 1;
						compList[i+1].rStart = compList[i].rEnd+1;
						tmpPos = compList[i+1].rStart + compList[i+1].seq.size() - 1;
						offset = tmpPos - compList[i+1].rEnd;
						compList[i+1].rEnd = tmpPos;
						updateCompPos(compList, i+2, offset);
						/////////////////
						// cerr<< "======== new compList ========" << endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].seq;
						// }
						// cerr<< endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].rStart << " " << compList[z].rEnd << endl;
						// }
						/////////////////
					}
				}
				else if(mapList_left.size()>0) // just left
				{
					// if(mapList_left.size()>1)
					// 	checkCnt++;
					// left
					if(mapList_left[0].rEnd >= compList[i].rEnd) // covers the whole un-corrected region
					{
						// cerr<< "case3" << endl;
						tmpStr = getQuerySubstr(mapList_left[0], compList[i].rStart, compList[i].rEnd);
						compList[i-1].seq += tmpStr + compList[i+1].seq;
						// adjust positions
						tmpPos = compList[i-1].rStart + compList[i-1].seq.size() - 1;
						offset = tmpPos - compList[i+1].rEnd;
						compList[i-1].rEnd = tmpPos;
						updateCompPos(compList, i+2, offset);
						compList.erase(compList.begin()+i, compList.begin()+i+2);
						i--;
						/////////////////
						// cerr<< "======== new compList ========" << endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].seq;
						// }
						// cerr<< endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].rStart << " " << compList[z].rEnd << endl;
						// }
						/////////////////
					}
					else // just partially covers the un-corrected region
					{
						// cerr<< "case4" << endl;
						tmpStr = getQuerySubstr(mapList_left[0], compList[i].rStart, mapList_left[0].rEnd);
						compList[i-1].seq += tmpStr;
						compList[i].seq = compList[i].seq.substr(mapList_left[0].rEnd-compList[i].rStart+1);
						// adjust positions
						compList[i-1].rEnd = compList[i-1].rStart + compList[i-1].seq.size() - 1;
						compList[i].rStart = compList[i-1].rEnd+1;
						tmpPos = compList[i].rStart + compList[i].seq.size() - 1;
						offset = tmpPos - compList[i].rEnd;
						compList[i].rEnd = tmpPos;
						updateCompPos(compList, i+1, offset);
						/////////////////
						// cerr<< "======== new compList ========" << endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].seq;
						// }
						// cerr<< endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].rStart << " " << compList[z].rEnd << endl;
						// }
						/////////////////
					}
				}
				else if(mapList_right.size()>0) // just right
				{
					// if(mapList_right.size()>1)
					// 	checkCnt++;
					// right
					if(mapList_right[0].rStart <= compList[i].rStart) // covers the whole un-corrected region
					{
						// cerr<< "case5" << endl;
						tmpStr = getQuerySubstr(mapList_right[0], compList[i].rStart, compList[i].rEnd);
						compList[i-1].seq += tmpStr + compList[i+1].seq;
						// adjust positions
						tmpPos = compList[i-1].rStart + compList[i-1].seq.size() - 1;
						offset = tmpPos - compList[i+1].rEnd;
						compList[i-1].rEnd = tmpPos;
						updateCompPos(compList, i+2, offset);
						compList.erase(compList.begin()+i, compList.begin()+i+2);
						i--;
						/////////////////
						// cerr<< "======== new compList ========" << endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].seq;
						// }
						// cerr<< endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].rStart << " " << compList[z].rEnd << endl;
						// }
						/////////////////
					}
					else
					{
						// cerr<< "case6" << endl;
						// right
						tmpStr = getQuerySubstr(mapList_right[0], mapList_right[0].rStart, compList[i].rEnd);
						compList[i+1].seq = tmpStr + compList[i+1].seq;
						compList[i].seq = compList[i].seq.substr(0, compList[i].seq.size()-(compList[i].rEnd-mapList_right[0].rStart+1));
						// adjust positions
						compList[i].rEnd = compList[i].rStart + compList[i].seq.size() - 1;
						compList[i+1].rStart = compList[i].rEnd+1;
						tmpPos = compList[i+1].rStart + compList[i+1].seq.size() - 1;
						offset = tmpPos - compList[i+1].rEnd;
						compList[i+1].rEnd = tmpPos;
						updateCompPos(compList, i+2, offset);
						/////////////////
						// cerr<< "======== new compList ========" << endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].seq;
						// }
						// cerr<< endl;
						// for(int z=0; z<compList.size(); z++)
						// {
						// 	cerr<< compList[z].rStart << " " << compList[z].rEnd << endl;
						// }
						/////////////////
					}
				}
				// exit(0);
			}
		}

		// *************************
		// build the corrected long read
		// *************************
		corrLR = "";
		for(i=0; i<compList.size(); i++)
		{
			corrLR += ( compList[i].isCorrected ? str2Upper(compList[i].seq) : str2Lower(compList[i].seq) );
		}

		// done with the correction!!!
		_cor_outBuffer[tid]<< ">" << currentLrId << " " << corrLR.size() << " " << tid << "\n"
			<< corrLR << "\n";
	}
	else
	{
		_cor_outBuffer[tid]<< ">" << currentLrId << " " << currentLrSeq.size() << " " << tid << "\n"
			<< currentLrSeq << "\n";
	}
}

void getComponents(string &lr, deque<component_t> &lst)
{
	lst.clear();
	int i=0;
	int s, e;
	while(i<lr.size())
	{
		s = i++;
		while(i<lr.size() && isupper(lr[i])==isupper(lr[s]))
			i++;
		e = i-1;
		lst.push_back(component_t(s, e, lr.substr(s, e-s+1), (isupper(lr[s])!=0 ? true : false)));
	}
}

void storeReadNIndex(deque<component_t> &compList, int tid)
{
	string com;
	com = "mkdir -p " + CURRENT_DIR + "/tmp_"+type2str<int>(tid);
	system(com.c_str());

	com = "rm -f " + CURRENT_DIR + "/tmp_"+type2str<int>(tid)+"/*";
	system(com.c_str());

	string lrSeq = "";
	for(int j=0; j<compList.size(); j++)
	{
		lrSeq += ( compList[j].isCorrected ? str2Upper(compList[j].seq) : str2Lower(compList[j].seq) );
	}
	ofstream foutLr((CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/lrSeq.fasta").c_str());
	if(foutLr.is_open()==false)
	{
		cerr<< "Could not open file: " << (CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/lrSeq.fasta") << endl;
		exit(0);
	}
	foutLr<< ">lrSeq" << "\n";
	foutLr<< lrSeq << endl;
	foutLr.close();

	// run BWA
	com = EXE_DIR+"/runBWAIndex.sh "+CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/lrSeq.fasta 1> /dev/null 2> /dev/null";
	system(com.c_str());
}

void extractOEAs(deque<sam_t> &i2pList, map<string, sam_t> &allMappings, set<string> &readList, int bPos, int ePos, bool isForward)
{
	int beginPos;
	int endPos;
	if(isForward == true)
	{
		beginPos = max(bPos, ePos-500);
		endPos = ePos;
	}
	else
	{
		beginPos = bPos;
		endPos = min(ePos, bPos+400);
	}

	readList.clear();
	for(int i=0; i<i2pList.size(); i++)
	{
		if(i2pList[i].rStart>=beginPos && i2pList[i].rStart<=endPos /*&& i2pList[i].score >= MIN_SCORE */&& i2pList[i].isForward==isForward)
		{
			string tmpStr = getMateId(i2pList[i].qName);
			if(allMappings.find(tmpStr)==allMappings.end())
			{
				// OEA read
				readList.insert(tmpStr);
			}
		}
		if(isForward==true && i2pList[i].rStart>=endPos-65 && i2pList[i].rStart<=endPos)
		{
			// edge read
			readList.insert(i2pList[i].qName);
		}
		if(isForward==false && i2pList[i].rStart>=beginPos && i2pList[i].rStart<=beginPos+5)
		{
			// edge read
			readList.insert(i2pList[i].qName);
		}
	}
}

void construtContigsFromOEA(set<string> &readList, deque<sam_t> &mapList, string fileName, map<string, string> &srList, int pivot, int tid)
{
	string com = "";
	// store oea short reads
	ofstream fout((CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/" + fileName).c_str());
	if(fout.is_open()==false)
	{
		cerr<< "Could not open file: " << (CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/" + fileName) << endl;
		exit(0);
	}
	for(set<string>::iterator it=readList.begin(); it!=readList.end(); it++)
	{
		fout<< ">" << *it << "\n";
		fout<< srList.find(*it)->second << "\n";
	}
	fout.close();

	// run Minia
	com  = EXE_DIR+"/runMinia.sh " + CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/ "+CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/" + fileName+" 43 1 " + CURRENT_DIR +"/tmp_"+type2str<int>(tid)+"/" + fileName + " 1> /dev/null 2> /dev/null";
	system(com.c_str());

	// run BWA
	com = EXE_DIR+"/runBWAMemPacbio.sh "+CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/"+fileName+".contigs.fa "+CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/lrSeq.fasta "+CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/" + fileName + ".sam 1> /dev/null 2> /dev/null";
	system(com.c_str());
	getAllSamEntries(CURRENT_DIR+"/tmp_"+type2str<int>(tid)+"/" + fileName + ".sam", mapList, pivot);
}

string getQuerySubstr(sam_t &a, int rs, int re)
{
	if(rs>re)
		return "";
    int cnt = 0;
    int qIndex = 0;
    int rIndex = a.rStart;
    int i;
	istringstream sin(a.cigar);
	int n;
	char c;
	string qStr="";
	while(sin>> n >> c)
	{
        switch (c)
        {
            case 'S':
            	if(cnt==0)
            		qIndex += n;
            	break;
            case 'M':
            case 'X':
            case '=':
            	for(i=0; i<n; i++)
            	{
            		if(rIndex>=rs && rIndex<=re)
            			qStr += a.read[qIndex];
            		qIndex++;
            		rIndex++;
            	}
            	break;
            case 'I':
            	for(i=0; i<n; i++)
            	{
            		if(rIndex>=rs && rIndex<=re)
            			qStr += a.read[qIndex];
            		qIndex++;
            	}
                break;
            case 'D':
            case 'N':
            case 'P':
            	rIndex += n;
                break;
        }
        cnt++;
    }
    return qStr;
}

void updateCompPos(deque<component_t> &compList, int index, int offset)
{
	int i;
	for(i=index; i<compList.size(); i++)
	{
		compList[i].rStart += offset;
		compList[i].rEnd += offset;
	}
}

string getMateId(string id)
{
	string mateId;
	int index = id.find_last_of('.');
	string prefix = id.substr(0, index+1);
	string suffix = id.substr(index+1);
	mateId = prefix + (suffix=="1" ? "2" : "1");
	return mateId;
}

void getAllSamEntries(string path, deque<sam_t> &lst, int pivot)
{
	ifstream finSam(path.c_str());
	if(finSam.is_open()==false)
	{
		cerr<< "Could not open file: " << path << endl;
		exit(0);
	}
	// get sam lines one by one
	string line, tmp;
	vector<string> field;
	lst.clear();
	int refLen, qStart, qEnd, qLen;
	while(getline(finSam, line))
	{
		if(line.empty() || line[0]=='@')
		{
			continue;
		}
		field.clear();
		istringstream sin(line);
		while(sin >> tmp)
		{
			field.push_back(tmp);
		}
		sam_t samTemp;
		samTemp.qName = field[0];
		samTemp.flag = str2type<int>(field[1]);
		samTemp.rName = field[2];
		samTemp.rStart = str2type<int>(field[3])-1;
		refLen = getMappedRefLen(field[5], qStart, qEnd, qLen);
		samTemp.rEnd = samTemp.rStart + refLen - 1;
		samTemp.qStart = qStart;
		samTemp.qEnd = qEnd;
		// samTemp.read = field[9].substr(qStart, qEnd-qStart+1);
		samTemp.read = field[9];
		samTemp.cigar = field[5];
		if(samTemp.rStart < pivot && samTemp.rEnd > pivot)
			lst.push_back(samTemp);
	}
}