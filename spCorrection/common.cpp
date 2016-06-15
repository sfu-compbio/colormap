/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#include "common.h"

// global in program
float 		MIN_IDENTITY = 60;
int 		THREAD_COUNT = 1;
int			THREAD_ID[255];
int 		MAX_READS_THREAD = 100;
int 		IS_DEBUG = 0;
string 		SR2LR_FILE = "";
string 		LR_FILE = "";
string 		OUTPUT_FILE = "";
string 		LOG_FILE = "";
int 		MIN_OVERLAP = 10-1; //minOverlap 10!!!

// global in common
char 		_com_revVal[128];
char 		_com_lowVal[128];
char 		_com_upVal[128];

void initCommon()
{
	_com_revVal['a'] = 't';
	_com_revVal['A'] = 'T';
	_com_revVal['c'] = 'g';
	_com_revVal['C'] = 'G';
	_com_revVal['g'] = 'c';
	_com_revVal['G'] = 'C';
	_com_revVal['t'] = 'a';
	_com_revVal['T'] = 'A';
	_com_revVal['n'] = 'n';
	_com_revVal['N'] = 'N';
	_com_revVal['-'] = '-';

	_com_lowVal['a'] = _com_lowVal['A'] = 'a';
	_com_lowVal['c'] = _com_lowVal['C'] = 'c';
	_com_lowVal['g'] = _com_lowVal['G'] = 'g';
	_com_lowVal['t'] = _com_lowVal['T'] = 't';
	_com_lowVal['n'] = _com_lowVal['N'] = 'n';

	_com_upVal['a'] = _com_upVal['A'] = 'A';
	_com_upVal['c'] = _com_upVal['C'] = 'C';
	_com_upVal['g'] = _com_upVal['G'] = 'G';
	_com_upVal['t'] = _com_upVal['T'] = 'T';
	_com_upVal['n'] = _com_upVal['N'] = 'N';

	int i;
	// cout<< "# threads: " << THREAD_COUNT << endl;
	for (i = 0; i < 255; i++)
		THREAD_ID[i] = i;
}

string str2Lower(string str)
{
	string lowerStr(str.size(), 'n');
	for(int i=0; i<str.size(); i++)
	{
		lowerStr[i] = _com_lowVal[str[i]];
	}
	return lowerStr;
}

string str2Upper(string str)
{
	string upperStr(str.size(), 'n');
	for(int i=0; i<str.size(); i++)
	{
		upperStr[i] = _com_upVal[str[i]];
	}
	return upperStr;
}

void revString(string &revStr, string &str)
{
	int n = str.size();
	revStr.assign(n, 'N');
	for(int i=0; i<n; i++)
		revStr[i] = str[n-i-1];
}

void revComplement(string &revStr, string &str)
{
	int n = str.size();
	revStr.assign(n, 'N');
	for(int i=0; i<n; i++)
		revStr[i] = _com_revVal[str[n-i-1]];
}

int parseCommandLine(int argc, char *argv[])
{
	int index, c;

	static struct option longOptions[] =
	{
		{"threads", 		required_argument, 		0, 				't'	},
		{"alignment", 		required_argument, 		0, 				'a'	},
		{"identity", 		required_argument, 		0, 				'i'	},
		{"longread", 		required_argument, 		0, 				'l'	},
		{"out", 			required_argument, 		0, 				'o'	},
		{"log", 			required_argument, 		0, 				'e'	},
		{"debug", 			no_argument, 			0, 				'd'	},
		{"help", 			no_argument, 			0, 				'h'	},
		{"version", 		no_argument, 			0, 				'v'	},
		{0,0,0,0}
	};

	while ( (c = getopt_long ( argc, argv, "t:a:i:l:o:e:dhv", longOptions, &index))!= -1 )
	{
		switch (c)
		{
			case 't':
				THREAD_COUNT = str2type<int>(optarg);
				if (THREAD_COUNT <= 0 || THREAD_COUNT > sysconf( _SC_NPROCESSORS_ONLN ))
					THREAD_COUNT = sysconf( _SC_NPROCESSORS_ONLN );
				break;
			case 'a':
				SR2LR_FILE = optarg;
				break;
			case 'i':
				MIN_IDENTITY = str2type<int>(optarg);
				break;
			case 'l':
				LR_FILE = optarg;
				break;
			case 'o':
				OUTPUT_FILE = optarg;
				break;
			case 'e':
				LOG_FILE = optarg;
				break;
			case 'd':
				IS_DEBUG = 1;
				break;
			case 'h':
				printHelp();
				return 0;
				break;
			case 'v':
				fprintf(stdout, "Version: %s\nBuild Date: %s\n", SPCORR_VERSION, BUILD_DATE);
				return 0;
				break;
			default:
				cerr << "[ERROR] run \"./spCorrection -h\" to get a list of acceptable arguments." << endl;
				return 0;
		}
	}

	if(LR_FILE == "")
	{
		cerr << endl;
		cerr << "[ERROR] option --longread/-l is required" << endl;
		cerr << " [NOTE] run \"./spCorrection -h\" to get a list of acceptable arguments." << endl;
		cerr << endl;
		return 0;
	}

	if(SR2LR_FILE == "")
	{
		cerr << endl;
		cerr << "[ERROR] option --alignment/-a is required" << endl;
		cerr << " [NOTE] run \"./spCorrection -h\" to get a list of acceptable arguments." << endl;
		cerr << endl;
		return 0;
	}

	if(MIN_IDENTITY<=0 || MIN_IDENTITY>100)
	{
		cerr << endl;
		cerr << "[ERROR] option --identity/-i requires an integer argument in range (0..100]" << endl;
		cerr << " [NOTE] run \"./spCorrection -h\" to get a list of acceptable arguments." << endl;
		cerr << endl;
		return 0;
	}

	return 1;
}

void printHelp()
{
	cerr << endl;
	cerr << "USAGE: spCorrection -l LR.fasta -a SR2LR.sam [options] > LR_corr.fasta" << endl;
	cerr << endl;
	cerr << "Required options:" << endl;
	cerr << "         -l STR        input set of long reads in fasta/q format" << endl;
	cerr << "         -a STR        mapping of short reads to long reads in sam format" << endl;
	cerr << "More options:" << endl;
	cerr << "         -t INT        number of threads [1]" << endl;
	cerr << "         -i INT        use top mappings with identity higher than INT [60]" << endl;
	cerr << "         -v            print version" << endl;
	cerr << "         -h            print this help" << endl;
	cerr << endl;
	cerr << "Note: In order to map short reads to long reads in m5 format use the following" << endl;
	cerr << "command: " << endl << endl;
	cerr << "      blasr SR.fasta LR.fasta -noSplitSubreads -m 5 -out sr2lr.m5 -nproc 8 -bestn 80 -nCandidates 80" << endl << endl;
	cerr << "then sort it based on the long reads ids:" << endl << endl;
	cerr << "      sort -k6,6V -o sr2lr.sorted.m5 sr2lr.m5" << endl;
	cerr << endl;
}
