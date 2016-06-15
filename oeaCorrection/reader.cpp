/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#include "reader.h"
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

gzFile 						_rd_finLr;
kseq_t 						*_rd_kseqLr;
ifstream 					_rd_finSr2Lr;
deque<read_t> 				_rd_lrList;
deque<read_t>::iterator 	_rd_lrItr;
pthread_mutex_t 			_rd_lrPosLock;
int 						_rd_maxReadChunk;
vector<string> 				_rd_lastField;

bool 						getNextSamEntry(ifstream &fin, vector<string> &fields);
void 						pushSamElem(deque<sam_t> &lst, map<string, sam_t> &mappingsAll, vector<string> &field, string longRead);

void initReader(string lrFile, string sr2LrFile)
{
	_rd_finLr = gzopen(lrFile.c_str(), "r");
	if(_rd_finLr==NULL)
	{
		cerr<< "Could not open file: " << lrFile << endl;
		exit(0);
	}
	_rd_kseqLr = kseq_init(_rd_finLr);

	_rd_finSr2Lr.open(sr2LrFile.c_str());
	if(_rd_finSr2Lr.is_open()==false)
	{
		cerr<< "Could not open file: " << sr2LrFile << endl;
		exit(0);
	}

	if(getNextSamEntry(_rd_finSr2Lr, _rd_lastField)==false) // read the first SAM entry
	{
		cerr<< "No sam entry in file: " << sr2LrFile << endl;
		exit(0);
	}

	_rd_maxReadChunk = THREAD_COUNT * MAX_READS_THREAD;
}

void finalizeReader()
{
	_rd_finSr2Lr.close();
	gzclose(_rd_finLr);
}

deque<read_t>* readChunk()
{
	_rd_lrList.clear();
	int cnt = 0;
	while(kseq_read(_rd_kseqLr) >= 0)
	{
		_rd_lrList.push_back(read_t(_rd_kseqLr->name.s, _rd_kseqLr->seq.s));
		cnt++;
		if(cnt == _rd_maxReadChunk)
			break;
	}
	if(_rd_lrList.size()>0)
	{
		_rd_lrItr = _rd_lrList.begin();
		return &_rd_lrList;
	}
	else
	{
		return NULL;
	}
}

read_t* getNextReadsInfo(deque<sam_t> &lstMap, map<string, sam_t>& allMap)
{
	read_t *p;
	lstMap.clear();
	pthread_mutex_lock(&_rd_lrPosLock);
	if(_rd_lrItr != _rd_lrList.end())
	{
		p = &(*_rd_lrItr);
		_rd_lrItr++;

		while(_rd_lastField[2]==p->id)
		{
			pushSamElem(lstMap, allMap, _rd_lastField, p->seq);
			if(getNextSamEntry(_rd_finSr2Lr, _rd_lastField)==false)
				break;
		}
	}
	else
	{
		p = NULL;
	}
	pthread_mutex_unlock(&_rd_lrPosLock);
	return p;
}

void getAllReads(string path, map<string, string> &readList)
{
	readList.clear();

	gzFile readFile;
	kseq_t *readObj;
	readFile = gzopen(path.c_str(), "r");
	if(readFile==NULL)
	{
		cerr<< "[getAllReads] ERROR: could not open file: " << path << endl;
		exit(0);
	}

	long long cnt = 0;
	readObj = kseq_init(readFile);
	while(kseq_read(readObj) >= 0)
	{
		cnt++;
		readList.insert(pair<string, string>(readObj->name.s, readObj->seq.s));
	}
	cerr << "[getAllReads] NOTE: " << cnt << " short reads loaded" << endl;
}

bool getNextSamEntry(ifstream &fin, vector<string> &fields)
{
	bool entryFound=false;
	string line;
	while(getline(fin, line))
	{
		if(line.empty()==false && line[0]!='@')
		{
			entryFound = true;
			break;
		}
		// cout<< "IGNORED!!!\t" << line << endl;
	}
	if(entryFound==false)
		return false;
	fields.clear();
	istringstream sin(line);
	string s;
	while(sin >> s)
	{
		fields.push_back(s);
	}
	return true;
}

void pushSamElem(deque<sam_t> &lst, map<string, sam_t> &mappingsAll, vector<string> &field, string longRead)
{
	int refLen, qStart, qEnd, qLen;
	sam_t mapTemp;
	mapTemp.qName = field[0];
	mapTemp.flag = str2type<int>(field[1]);
	mapTemp.isForward = ((mapTemp.flag & 16)==0 ? 1 : 0);
	mapTemp.rName = field[2];
	mapTemp.rStart = str2type<int>(field[3])-1;
	refLen = getMappedRefLen(field[5], qStart, qEnd, qLen);
	mapTemp.rEnd = mapTemp.rStart + refLen - 1;
	mapTemp.qStart = qStart;
	mapTemp.qEnd = qEnd;
	mapTemp.read = field[9];
	mapTemp.cigar = field[5];
	mapTemp.numMatch = getAlignmentMatchNum(mapTemp.cigar, mapTemp.read, longRead, mapTemp.rStart);
	mapTemp.identity = mapTemp.numMatch/(float)(mapTemp.qEnd-mapTemp.qStart+1);
	lst.push_back(mapTemp);
	mappingsAll.insert(pair<string, sam_t>(field[0], mapTemp));
}