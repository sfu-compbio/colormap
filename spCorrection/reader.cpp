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

bool 						getNextSamEntry(ifstream &fin, vector<string> &field);
void 						pushSamElem(deque<sam_t> &lst, vector<string> &field, string longRead);
int 						getMappedRefLen(string cigar, int &qStart, int &qEnd, int &qLen);
float 						getAlignmentIdentity(string cigar, string qSeq, string rSeq, int rStart);

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
	if(getNextSamEntry(_rd_finSr2Lr, _rd_lastField)==false) // read the first sam entry
	{
		cerr<< "No sam entries in file: " << sr2LrFile << endl;
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
		_rd_lrList.push_back(read_t(_rd_kseqLr->name.s, str2Lower(_rd_kseqLr->seq.s)));
		cnt++;
		if(cnt == _rd_maxReadChunk)
			break;
	}
	// cout<< cnt << endl;
	// char dummy;
	// cin >> dummy;
	// exit(0);
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

read_t* getNextReadsInfo(deque<sam_t> &lst)
{
	read_t *p;
	lst.clear();
	pthread_mutex_lock(&_rd_lrPosLock);
	if(_rd_lrItr != _rd_lrList.end())
	{
		p = &(*_rd_lrItr);
		_rd_lrItr++;

		while(_rd_lastField[2]==p->id)
		{
			pushSamElem(lst, _rd_lastField, p->seq);
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

bool getNextSamEntry(ifstream &fin, vector<string> &field)
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
	}
	if(entryFound==false)
		return false;
	field.clear();
	istringstream sin(line);
	string s;
	while(sin >> s)
	{
		field.push_back(s);
	}
	return true;
}

void pushSamElem(deque<sam_t> &lst, vector<string> &field, string longRead)
{
	int refLen, qStart, qEnd, qLen;
	sam_t mapTemp;
	mapTemp.qName = field[0];
	mapTemp.flag = str2type<int>(field[1]);
	mapTemp.rName = field[2];
	mapTemp.rStart = str2type<int>(field[3])-1;
	refLen = getMappedRefLen(field[5], qStart, qEnd, qLen);
	mapTemp.rEnd = mapTemp.rStart + refLen - 1;
	mapTemp.qStart = qStart;
	mapTemp.qEnd = qEnd;
	mapTemp.read = str2Lower(field[9]);
	mapTemp.quality = field[10];
	mapTemp.cigar = field[5];
	mapTemp.identity = getAlignmentIdentity(mapTemp.cigar, mapTemp.read, longRead, mapTemp.rStart);
	if(mapTemp.identity >= MIN_IDENTITY)
		lst.push_back(mapTemp);
}

deque<read_t> getAllReads(string path)
{
	deque<read_t> readList;

	gzFile readFile;
	kseq_t *readObj;
	readFile = gzopen(path.c_str(), "r");
	if(readFile==NULL)
	{
		cerr<< "Could not open file: " << path << endl;
		exit(0);
	}

	readObj = kseq_init(readFile);
	while(kseq_read(readObj) >= 0)
	{
		readList.push_back(read_t(readObj->name.s, readObj->seq.s));
	}

	return readList;
}

int getMappedRefLen(string cigar, int &qStart, int &qEnd, int &qLen)
{
    int contLen = 0;
    int readLen = 0;
    int beginClip=0, endClip=0;
    int cnt = 0;
	istringstream sin(cigar);
	int n;
	char c;
	while(sin>> n >> c)
	{
        switch (c)
        {
            case 'S':
            case 'H':
            	if(cnt==0)
            		beginClip = n;
            	break;
            case 'M':
            case 'X':
            case '=':
                contLen += n;
                readLen += n;
                break;
            case 'I':
                readLen += n;
                break;
            case 'D':
            case 'N':
            case 'P':
                contLen += n;
                break;
        }
        cnt++;
    }
    if(c=='S' || c=='H')
    	endClip = n;
    qStart = beginClip;
    qEnd = beginClip+readLen-1;
    qLen = beginClip+readLen+endClip;
    return contLen;
}

float getAlignmentIdentity(string cigar, string qSeq, string rSeq, int rStart)
{
    int cnt = 0;
    int nMatch = 0;
    int qIndex = 0;
    int rIndex = rStart;
    int i;
	istringstream sin(cigar);
	int n;
	char c;
	while(sin>> n >> c)
	{
        switch (c)
        {
            case 'S':
            	if(cnt==0)
            	{
            		qIndex += n;
            	}
            	break;
            // case 'H':
            case 'M':
            	for(i=0; i<n; i++)
            	{
            		if(qSeq[qIndex]==rSeq[rIndex])
            		{
            			nMatch++;
            		}
            		qIndex++;
            		rIndex++;
            	}
            	break;
            case 'X':
            	for(i=0; i<n; i++)
            	{
            		qIndex++;
            		rIndex++;
            	}
            	break;
            case '=':
            	for(i=0; i<n; i++)
            	{
        			nMatch++;
            		qIndex++;
            		rIndex++;
            	}
                break;
            case 'I':
            	qIndex += n;
                break;
            case 'D':
            case 'N':
            case 'P':
            	rIndex += n;
                break;
        }
        cnt++;
    }
    return (float)nMatch/qSeq.size()*100;
}