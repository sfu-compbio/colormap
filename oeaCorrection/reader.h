/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#ifndef __READER__
#define __READER__

#include "common.h"

class read_t
{
public:
	string id;
	string seq;
	//
	read_t(string readId, string readSeq)
	{
		id = readId;
		seq = readSeq;
	}
};

typedef struct
{
	string qName;
	int flag;
	bool isForward;
	string rName;
	int rStart;
	int rEnd;
	int qStart;
	int qEnd;
	int numMatch;
	float identity;
	string read;
	string cigar;
} sam_t;

void 				initReader(string lrFile, string sr2LrFile);
void 				finalizeReader();
deque<read_t>* 		readChunk();
read_t* 			getNextReadsInfo(deque<sam_t> &lstMap, map<string, sam_t>& allMap);
void 				getAllReads(string path, map<string, string> &readList);

#endif //__READER__