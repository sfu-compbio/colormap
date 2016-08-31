/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#ifndef __CORRECT__
#define __CORRECT__

#include "common.h"

class component_t
{
public:
	int rStart;
	int rEnd;
	string seq;
	bool isCorrected;
	//
	component_t(int s, int e, string str, bool isCorr)
	{
		rStart = s;
		rEnd = e;
		seq = str;
		isCorrected = isCorr;
	}
};

class edge_t
{
public:
	int source;
	int dest;
	// string srcStr;
	// string destStr;
	double weight; // edit distance of the non-overlapped part of the second read
	// string suffix;
	int suffixPos;
	int srcPos;
	int destPos;
	//
	edge_t()
	{
		source = -1;
		dest = -1;
		weight = -1;
	}
};

void initCorrect();
void finalizeCorrect();
void correctLrMT();

#endif //__CORRECT__