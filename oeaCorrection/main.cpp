/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#include "common.h"
#include "reader.h"
#include "correct.h"
#include "kseq.h"

using namespace std;

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[])
{
	//
	if (!parseCommandLine(argc, argv))
	{
		return 1;
	}

	initCommon();
	initReader(LR_FILE, SR2LR_FILE);
	initCorrect();

	deque<read_t> *lrListP;
	while((lrListP=readChunk())!=NULL)
	{
		correctLrMT();
	}

	finalizeCorrect();
	finalizeReader();

	return 0;
}