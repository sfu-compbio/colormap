/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#include "correct.h"
#include "reader.h"
#include <pthread.h>

pthread_t					*_cor_threads = NULL;
pthread_mutex_t 			_cor_printLock;
ostringstream 				*_cor_outBuffer;
const double 				_cor_max_weight = std::numeric_limits<double>::infinity();
int 						MAX_MISMATCH = 1;

void* 						correctLrST(void *idp);
void 						correctSingleLr(int tid, read_t *rp, deque<sam_t> &i2pList);
bool 						compareSam(const sam_t &a, const sam_t &b);
string 						processConnectedComponent(deque<sam_t> &nodeList, string rSeq, int &startInRef, int &endInRef);
float 						calcEdit(sam_t &mapObj, string rSeq, int startPos);
void 						dijkstraDistance(deque<deque<edge_t> > &G, int source, vector<double> &min_distance, /*vector<int> &previous, */vector<edge_t> &prev_edge);
deque<int> 					dijkstraGetPath(int vertex, const vector<int> &previous);
deque<edge_t> 				dijkstraGetPath2(int vertex, const vector<edge_t> &previous);
int 						readSuffixPos(sam_t &aln1, int qs, sam_t &aln2, int rs);
int 						readPrefixPos(sam_t &a, int rs);
int 						findSuffixPos(sam_t &mapObj2, int suffixRefPos, float &weight);
bool 						checkOverlap(sam_t &mapObj1, sam_t &mapObj2, int qPos, string &str1, string &str2);


void initCorrect()
{
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
	// cout<< "in correctLrMT" << endl;
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
	deque<sam_t> m5List;
	while((rp=getNextReadsInfo(m5List)) != NULL)
	{
		// pthread_mutex_lock(&_rd_printLock);
		// cout << rp->id << "\t" << m5List.size() << "\t" << tid << endl;
		// pthread_mutex_unlock(&_rd_printLock);

		correctSingleLr(tid, rp, m5List);
	}
}

void correctSingleLr(int tid, read_t *rp, deque<sam_t> &i2pList)
{
	deque<sam_t> nodeList;
	string currentLrId = rp->id;
	string &currentLrSeq = rp->seq;
	// string unCorrLR = str2Lower(currentLrSeq);
	string unCorrLR = currentLrSeq;

	// cerr<< ">>>" << currentLrId << endl << endl << endl;
	// cerr<< ">>>" << currentLrId << endl;

	if(i2pList.size() > 0)
	{
		sort(i2pList.begin(), i2pList.end(), compareSam);

		// cerr<< ">>>" << currentLrId << " " << i2pList.size() << endl;
		// for(int i=0; i<i2pList.size(); i++)
		// 	cerr<< "\t" << i2pList[i].qName << "\t" <<  i2pList[i].rName << "\t" << i2pList[i].rStart << "\t"
		// 	<< i2pList[i].rEnd << "\t" << i2pList[i].qStart << "\t" << i2pList[i].qEnd << "\t" << i2pList[i].identity << endl;

		int i;
		int sPos;
		int ePos;
		int startInRef, endInRef;
		string assembledString;
		string corrLR = "";
		int lastPosInLR = 0;
		bool firstComponent = true;
		for(i=0; i<i2pList.size(); i++) // find connected components
		{
			if(firstComponent)
			{
				nodeList.push_back(i2pList[i]);
				sPos = i2pList[i].rStart;
				ePos = i2pList[i].rEnd;
				firstComponent = false;
			}
			else
			{
				if(i2pList[i].rStart <= ePos - MIN_OVERLAP) // overlap
				{
					nodeList.push_back(i2pList[i]);
					ePos = i2pList[i].rEnd;
				}
				else // no overlap with minOverlap
				{
					// process the current connected component
					assembledString = processConnectedComponent(nodeList, currentLrSeq, startInRef, endInRef);
					// assembledString = callConsensus(nodeList, startInRef, endInRef);
					if(assembledString != "")
					{
						if(startInRef<lastPosInLR)
						{
							if(endInRef>=lastPosInLR && lastPosInLR-startInRef<assembledString.size())
							{
								corrLR += assembledString.substr(lastPosInLR-startInRef);
								lastPosInLR = endInRef+1;
							}
						}
						else
						{
							corrLR += unCorrLR.substr(lastPosInLR, startInRef-lastPosInLR) + assembledString;
							lastPosInLR = endInRef+1;
						}
					}
					// create the new connected component
					nodeList.clear();
					nodeList.push_back(i2pList[i]);
					sPos = i2pList[i].rStart;
					ePos = i2pList[i].rEnd;
				}
			}
		}
		// process the last connected component
		assembledString = processConnectedComponent(nodeList, currentLrSeq, startInRef, endInRef);
		// assembledString = callConsensus(nodeList, startInRef, endInRef);
		if(assembledString != "")
		{
			if(startInRef<lastPosInLR)
			{
				if(endInRef>=lastPosInLR && lastPosInLR-startInRef<assembledString.size())
				{
					corrLR += assembledString.substr(lastPosInLR-startInRef);
					lastPosInLR = endInRef+1;
				}
			}
			else
			{
				corrLR += unCorrLR.substr(lastPosInLR, startInRef-lastPosInLR) + assembledString;
				lastPosInLR = endInRef+1;
			}
		}
		if(lastPosInLR < unCorrLR.size())
			corrLR += unCorrLR.substr(lastPosInLR);

		// done with the correction!!!
		// cout<< ">" << currentLrId << " " << unCorrLR.size() << "\n"
		// 	<< unCorrLR << "\n";
		_cor_outBuffer[tid]<< ">" << currentLrId << /*"_corr" <<*/ " " << corrLR.size() << "\n"
			<< corrLR << "\n";
		// cout << "\n";
	}
	else
	{
		_cor_outBuffer[tid]<< ">" << currentLrId << /*"_corr" <<*/ " " << unCorrLR.size() << "\n"
			<< unCorrLR << "\n";
	}
}

bool compareSam(const sam_t &a, const sam_t &b)
{
	if(a.rStart < b.rStart)
		return true;
	else if(a.rStart > b.rStart)
		return false;
	else if(a.rEnd < b.rEnd)
		return true;
	else
		return false;
}

string processConnectedComponent(deque<sam_t> &nodeList, string rSeq, int &startInRef, int &endInRef)
{
	int i, j, k, z;
	int p1;
	float w1;
	int offset[10] = {0, 1, -1, 2, -2, 3, -3, 4, -4};
	int startNode, endNode;
	int endNodePos = -1;

	// construct the graph now
	deque<deque<edge_t> > graph(nodeList.size());
	int numOfEdges = 0;
	for(j=0; j<nodeList.size()-1; j++)
	{
		for(k=j+1; k<nodeList.size(); k++)
		{
			// if(nodeList[k].rStart <= nodeList[j].rEnd - MIN_OVERLAP && nodeList[j].rEnd < nodeList[k].rEnd) // add an edge just in case of overlap (with at least MIN_OVERLAP bases)
			if(nodeList[k].rStart + MIN_OVERLAP <= nodeList[j].rEnd && nodeList[j].rEnd < nodeList[k].rEnd && nodeList[j].rStart < nodeList[k].rStart) // add an edge just in case of overlap (with at least MIN_OVERLAP bases)
			{
				edge_t e1;
				string tmpStr1, tmpStr2;
				e1.dest = k;
				p1 = readPrefixPos(nodeList[k], nodeList[j].rEnd);
				w1 = calcEdit(nodeList[k], rSeq, nodeList[j].rEnd);
				e1.source = j;

				for(z=0; z<7; z++)
				{
					if(checkOverlap(nodeList[j], nodeList[k], p1+offset[z], tmpStr1, tmpStr2))
					{
						if(p1+offset[z] >= nodeList[k].qEnd)
							continue;
						// e1.weight = w1+abs(offset[z]);
						e1.weight = w1;
						// e1.weight = (w1+abs(offset[z]))/(nodeList[k].read.size()-(p1+offset[z])+1);
						
						// e1.suffix = nodeList[k].read.substr(p1+offset[z], nodeList[k].qEnd-(p1+offset[z])+1);
						// e1.suffix = tmpStr2.substr(p1+offset[z], nodeList[k].qEnd-(p1+offset[z])+1);
						e1.suffixPos = p1+offset[z];
						
						e1.srcStr = tmpStr1;
						e1.destStr = tmpStr2;

						graph[j].push_back(e1);
						numOfEdges++;
						if(numOfEdges==1)
						{
							startNode = j;
						}
						if(nodeList[k].rEnd > endNodePos)
							endNode = k;
						// cerr<< "+" << endl;
						break;
					}
				}
				// if(z==7)
				// {
				// 	cerr<< "*" << endl;
				// }
			}
		}
	}

	if(numOfEdges==0)
	{
		// ERROR!!!
		// cerr<< "@" << endl;
		return "";
	}

	// cerr<< "+ component "<< nodeList[startNode].rStart << " " << nodeList[endNode].rEnd << endl;

	// int startNode = 0;
	// int endNode = nodeList.size()-1;
	vector<double> min_distance;
	vector<int> previous;
	vector<edge_t> prev_edge;
	dijkstraDistance(graph, startNode, min_distance, /*previous*/prev_edge);
	string spelledStr = "";
	// int tmpPos1, tmpPos2;
	if(min_distance[endNode] != _cor_max_weight)
	{
		deque<edge_t> path = dijkstraGetPath2(endNode, prev_edge);

		//
		// int indent = 0;
		// int indentSoFar = 0;
		// for(j=1; j<path.size(); j++)
		// {
		// 	cerr << string(indentSoFar, ' ') << path[j].srcStr << endl;
		// 	cerr << string(indentSoFar, ' ') << nodeList[path[j].source].quality << endl;
		// 	indent = path[j].srcStr.size()-path[j].suffixPos;
		// 	cerr << string(indentSoFar, ' ') << string(indent, ' ') << path[j].destStr << endl;
		// 	cerr << string(indentSoFar, ' ') << string(indent, ' ') << nodeList[path[j].dest].quality << endl << endl;
		// 	indentSoFar += indent;
		// }
		//

		startInRef = nodeList[path[1].source].rStart;
		// spelledStr += nodeList[path[1].source].read;
		spelledStr += path[1].srcStr;

		for(j=1; j<path.size()-1; j++)
		{
			// spelledStr += path[j].suffix;
			spelledStr += path[j+1].srcStr.substr(path[j].suffixPos, nodeList[path[j+1].source].qEnd - path[j].suffixPos + 1);
		}
		spelledStr += path[j].destStr.substr(path[j].suffixPos, nodeList[path[j].dest].qEnd - path[j].suffixPos + 1);
		endInRef = nodeList[path[path.size()-1].dest].rEnd;
		return str2Upper(spelledStr);
	}
	else
	{
		// ERROR!!!
		// cout<< "ERROR!!!" << endl;
		// exit(0);
		return "";
	}
}

float calcEdit(sam_t &mapObj, string rSeq, int startPos)
{
    int cnt = 0;
    int qIndex = 0;
    int rIndex = mapObj.rStart;
    string qSeq = mapObj.read;
	int edit = 0;
	int remLen = 0;
    int i;
	istringstream sin(mapObj.cigar);
	int n;
	char c;
	while(sin>> n >> c)
	{
        switch (c)
        {
            case 'S':
            	if(cnt==0)
            	{
	            	for(i=0; i<n; i++)
	            	{
	            		if(rIndex > startPos)
	            			remLen++;
	            		qIndex++;
	            	}
            	}
            	break;
            // case 'H':
            case 'M':
            	for(i=0; i<n; i++)
            	{
            		if(rIndex > startPos)
            		{
	            		if(tolower(qSeq[qIndex])!=tolower(rSeq[rIndex]))
	            			edit++;
	            		remLen++;
            		}
            		qIndex++;
            		rIndex++;
            	}
            	break;
            case 'X':
            	for(i=0; i<n; i++)
            	{
            		if(rIndex > startPos)
            		{
            			edit++;
            			remLen++;
            		}
            		qIndex++;
            		rIndex++;
            	}
            	break;
            case '=':
            	for(i=0; i<n; i++)
            	{
            		if(rIndex > startPos)
            		{
            			remLen++;
            		}
            		qIndex++;
            		rIndex++;
            	}
                break;
            case 'I':
            	if(rIndex > startPos)
            	{
            		edit += n;
            		remLen += n;
            	}
            	qIndex += n;
                break;
            case 'D':
            case 'N':
            case 'P':
            	for(i=0; i<n; i++)
            	{
            		if(rIndex > startPos)
            			edit++;
            		rIndex++;
            	}
                break;
        }
        cnt++;
    }
	// cout<< "\t+\t" << startPos << "\t" << edit << "\t" << remLen << endl;
    return (float)edit/remLen;
}

void dijkstraDistance(deque<deque<edge_t> > &G, int source, vector<double> &min_distance, /*vector<int> &previous, */vector<edge_t> &prev_edge)
{
	int n = G.size();
	min_distance.clear();
	min_distance.resize(n, _cor_max_weight);
	min_distance[source] = 0;
	// previous.clear();
	// previous.resize(n, -1);
	prev_edge.clear();
	prev_edge.resize(n);
	set<pair<double, int> > vertex_queue;
	vertex_queue.insert(make_pair(min_distance[source], source));

	while (!vertex_queue.empty()) 
	{
		double dist = vertex_queue.begin()->first;
		int u = vertex_queue.begin()->second;
		vertex_queue.erase(vertex_queue.begin());

		// Visit each edge exiting u
		const deque<edge_t> &neighbors = G[u];
		for (deque<edge_t>::const_iterator neighbor_iter = neighbors.begin(); neighbor_iter != neighbors.end(); neighbor_iter++)
		{
			int v = neighbor_iter->dest;
			double weight = neighbor_iter->weight;
			double distance_through_u = dist + weight;
			if (distance_through_u < min_distance[v])
			{
				vertex_queue.erase(make_pair(min_distance[v], v));

				min_distance[v] = distance_through_u;
				// previous[v] = u;
				prev_edge[v] = *neighbor_iter;
				vertex_queue.insert(make_pair(min_distance[v], v));
			}
		}
	}
}

deque<int> dijkstraGetPath(int vertex, const vector<int> &previous)
{
	deque<int> path;
	for ( ; vertex != -1; vertex = previous[vertex])
		path.push_front(vertex);
	return path;
}

/*
aln1  ==================
aln2               ===================
*/

int readSuffixPos(sam_t &aln1, int qs, sam_t &aln2, int rs)
{
	// cerr<< aln1.read << endl;
	// cerr<< string(qs, ' ') << aln2.read << endl;

	int z;
	int len = aln1.read.size()-qs;
	int err=0;
	for(z=0; z<len; z++)
	{
		if(aln1.read[z+qs]!=aln2.read[z])
			err++;
	}
	if(err==0)
	{
		// cerr<< z << endl;
		return z;
	}
	// cerr<< err << endl << "=================" << endl;

	if(rs<aln2.rStart)
	{
		cerr<< "Error in getReadSuffix(): Wrong start on the reference!!!" << endl;
		exit(0);
	}

	if(rs == aln2.rStart)
	{
		return 1;
	}

    int cnt = 0;
    int qIndex = 0;
    int rIndex = aln2.rStart;
    int i;
	istringstream sin(aln2.cigar);
	int n;
	char c;
	while(sin>> n >> c)
	{
        switch (c)
        {
            case 'S':
            	if(cnt==0)
            		qIndex += n;
            	break;
            // case 'H':
            case 'M':
            case 'X':
            case '=':
            	for(i=0; i<n; i++)
            	{
            		qIndex++;
            		rIndex++;
            		if(rIndex == rs)
						return qIndex+1;
            	}
            	break;
            case 'I':
            	qIndex += n;
                break;
            case 'D':
            case 'N':
            case 'P':
            	for(i=0; i<n; i++)
            	{
            		rIndex++;
            		if(rIndex == rs)
						return qIndex+1;
            	}
                break;
        }
        cnt++;
    }
    return qIndex+1;
}

int readPrefixPos(sam_t &aln, int rs)
{
	if(rs<aln.rStart)
	{
		cerr<< "Error in getReadSuffix(): Wrong start on the reference!!!" << endl;
		exit(0);
	}

	if(rs == aln.rStart)
	{
		return 0;
	}

    int cnt = 0;
    int qIndex = 0;
    int rIndex = aln.rStart;
    int i;
	istringstream sin(aln.cigar);
	int n;
	char c;
	while(sin>> n >> c)
	{
        switch (c)
        {
            case 'S':
            	if(cnt==0)
            		qIndex += n;
            	break;
            // case 'H':
            case 'M':
            case 'X':
            case '=':
            	for(i=0; i<n; i++)
            	{
            		qIndex++;
            		rIndex++;
            		if(rIndex == rs)
						return qIndex;
            	}
            	break;
            case 'I':
            	qIndex += n;
                break;
            case 'D':
            case 'N':
            case 'P':
            	for(i=0; i<n; i++)
            	{
            		rIndex++;
            		if(rIndex == rs)
						return qIndex;
            	}
                break;
        }
        cnt++;
    }
    return qIndex;
}

// int findSuffixPos(sam_t &mapObj2, int suffixRefPos, float &weight)
// {
// 	int rPos = mapObj2.rEnd;
// 	int qPos = mapObj2.read.size()-1;
// 	int transLen = mapObj2.qTrans.size();
// 	int edit = 0;
// 	int remLen = 0;
// 	for(int i=transLen-1; i>=0 && rPos>suffixRefPos; i--)
// 	{
// 		if(mapObj2.mTrans[i]=='|')
// 		{
// 			rPos--;
// 			qPos--;
// 		}
// 		else // mapObj.mTrans[i]=='*'
// 		{
// 			if(mapObj2.qTrans[i]=='-') // deletion
// 			{
// 				edit++;
// 				rPos--;
// 			}
// 			else if(mapObj2.rTrans[i]=='-') // insertion
// 			{
// 				edit++;
// 				qPos--;
// 			}
// 			else // mis-match
// 			{
// 				edit++;
// 				rPos--;
// 				qPos--;
// 			}
// 		}
// 	}
// 	weight = edit;
// 	return qPos+1;
// }

bool checkOverlap(sam_t &mapObj1, sam_t &mapObj2, int qPos, string &str1, string &str2)
{
	if(qPos >= mapObj2.read.size())
		return false;
	str1 = mapObj1.read;
	str2 = mapObj2.read;
	int l = mapObj1.qEnd;
	int cnt = 0;
	// cerr<< "+++" << endl;
	// cerr<< mapObj1.qName << " " << mapObj1.qStart << " " << mapObj1.qEnd << " " << mapObj1.rStart << " " << mapObj1.rEnd << endl;
	// cerr<< mapObj2.qName << " " << mapObj2.qStart << " " << mapObj2.qEnd << " " << mapObj2.rStart << " " << mapObj2.rEnd << endl;
	for(int p = qPos-1; l>=mapObj1.qStart && p>=mapObj2.qStart; p--, l--)
	{
		// cerr<< l << "/" << mapObj1.read.size()-1 << "\t" << p << "/" << mapObj2.read.size()-1 << endl;
		if(mapObj1.read[l]!=mapObj2.read[p])
		{
			if(mapObj1.quality[l] > mapObj2.quality[p])
				str2[p] = str1[l];
			else if(mapObj1.quality[l] < mapObj2.quality[p])
				str1[l] = str2[p];
			cnt++;
		}
	}

    return (cnt<=1);
    // return (cnt<=2);
    // return (cnt<1);
}

deque<edge_t> dijkstraGetPath2(int vertex, const vector<edge_t> &previous)
{
	deque<edge_t> path;
	for ( ; vertex != -1; vertex = previous[vertex].source)
		path.push_front(previous[vertex]);
	return path;
}