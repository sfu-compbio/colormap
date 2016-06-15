/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#ifndef __COMMON__
#define __COMMON__

#include <stdlib.h>
#include <string.h>
#include <string>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>
#include <math.h>
#include <stdint.h>
#include <limits>
#include <zlib.h>
#include <unistd.h>
#include <getopt.h>

using namespace std;

extern float 		MIN_IDENTITY;
extern int 			THREAD_COUNT;
extern int			THREAD_ID[255];
extern int 			MAX_READS_THREAD;
extern int 			EXTEND_CORR;
extern int 			IS_DEBUG;
extern string 		SR2LR_FILE;
extern string 		LR_FILE;
extern string 		SR_FILE;
extern string 		OUTPUT_FILE;
extern string 		LOG_FILE;
extern string 		CURRENT_DIR;
extern string 		EXE_DIR;

void 				initCommon();
string 				str2Lower(string str);
string 				str2Upper(string str);
void 				revString(string &str, string &revStr);
void 				revComplement(string &str, string &revStr);
int 				parseCommandLine (int argc, char *argv[]);
void 				printHelp();
int 				getMappedRefLen(string &cigar, int &qStart, int &qEnd, int &qLen);
int 				getAlignmentMatchNum(string cigar, string qSeq, string rSeq, int rStart);

template <typename T>
T str2type(string str)
{
	T n;
	istringstream sin(str);
	sin >> n;
	return n;
}

template <typename T>
string type2str(T v)
{
	ostringstream sout;
	sout << v;
	return sout.str();
}

#endif //__COMMON__