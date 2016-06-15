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
	int flag;
	string qName;
	int qStart;
	int qEnd;
	string rName;
	int rStart;
	int rEnd;
	float identity;
	string read;
	string quality;
	string cigar;
} sam_t;

void 				initReader(string lrFile, string sr2LrFile);
void 				finalizeReader();
deque<read_t>* 		readChunk();
read_t*				getNextReadsInfo(deque<sam_t> &lst);
deque<read_t> 		getAllReads(string path);

#endif //__READER__