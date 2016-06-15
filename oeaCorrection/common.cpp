/*
 * Author: Ehsan Haghshenas (ehaghshe[at]sfu[dot]ca)
 */

#include "common.h"

// global in program
float 		MIN_IDENTITY = 75;
int 		THREAD_COUNT = 1;
int			THREAD_ID[255];
int 		MAX_READS_THREAD = 20;
int 		EXTEND_CORR = 0;
int 		IS_DEBUG = 0;
string 		SR2LR_FILE = "";
string 		LR_FILE = "";
string 		SR_FILE = "";
string 		OUTPUT_FILE = "";
string 		LOG_FILE = "";
string 		CURRENT_DIR = "";
string 		EXE_DIR = "";

// global in common
char 		_com_revVal[128];
char 		_com_lowVal[128];
char 		_com_upVal[128];

string 		getCurrentDirectory();
string 		getExeDirectory();

void initCommon()
{
	CURRENT_DIR = getCurrentDirectory();
	EXE_DIR = getExeDirectory();

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
		{"longRead", 		required_argument, 		0, 				'l'	},
		{"shortRead", 		required_argument, 		0, 				's'	},
		{"out", 			required_argument, 		0, 				'o'	},
		{"log", 			required_argument, 		0, 				'e'	},
		{"debug", 			no_argument, 			0, 				'd'	},
		{"help", 			no_argument, 			0, 				'h'	},
		{"version", 		no_argument, 			0, 				'v'	},
		{0,0,0,0}
	};

	while ( (c = getopt_long ( argc, argv, "t:a:l:s:o:e:dhv", longOptions, &index))!= -1 )
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
			case 'l':
				LR_FILE = optarg;
				break;
			case 's':
				SR_FILE = optarg;
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
				fprintf(stdout, "Version: %s\nBuild Date: %s\n", OEACORR_VERSION, BUILD_DATE);
				return 0;
				break;
			default:
				cerr << "[ERROR] run \"./oeaCorrection -h\" to get a list of acceptable arguments." << endl;
				return 0;
		}
	}

	if(SR_FILE == "")
	{
		cerr << endl;
		cerr<< "[ERROR] option --shortRead/-s is required" << endl;
		cerr << " [NOTE] run \"./oeaCorrection -h\" to get a list of acceptable arguments." << endl;
		cerr << endl;
		return 0;
	}

	if(LR_FILE == "")
	{
		cerr << endl;
		cerr<< "[ERROR] option --longRead/-l is required" << endl;
		cerr << " [NOTE] run \"./oeaCorrection -h\" to get a list of acceptable arguments." << endl;
		cerr << endl;
		return 0;
	}

	if(SR2LR_FILE == "")
	{
		cerr << endl;
		cerr<< "[ERROR] option --alignment/-a is required" << endl;
		cerr << " [NOTE] run \"./oeaCorrection -h\" to get a list of acceptable arguments." << endl;
		cerr << endl;
		return 0;
	}

	return 1;
}

void printHelp()
{
	cerr << endl;
	cerr << "USAGE: oeaCorrection -s SR.fasta -l LR_corr.fasta -a SR2LR.sam [options] > LR_oea.fasta" << endl;
	cerr << endl;
	cerr << "Required options:" << endl;
	cerr << "         -s STR        input set of short reads in fasta/q format" << endl;
	cerr << "         -l STR        input set of long reads corrected by spCorrection" << endl;
	cerr << "         -a STR        mapping of short reads to long reads in m5 format" << endl;
	cerr << "More options:" << endl;
	cerr << "         -t INT        number of threads [1]" << endl;
	cerr << "         -v            print version" << endl;
	cerr << "         -h            print this help" << endl;
	cerr << endl;
	cerr << "Note: In order to map short reads to long reads (corrected by spCorrection)" << endl;
	cerr << "use the following command: " << endl << endl;
	cerr << "      utils/bwa-colormap/bwa index LR_corr.fasta" << endl;
	cerr << "      utils/bwa-colormap/bwa mem -t 8 -aY LR_corr.fasta illumina.fastq > sr2corr.sam" << endl << endl;
	cerr << "then sort it based on the long reads ids:" << endl << endl;
	cerr << "      utils/sortSam.sh sr2corr.sam sr2corr.sorted.sam" << endl;
	cerr << endl;
}

int getMappedRefLen(string &cigar, int &qStart, int &qEnd, int &qLen)
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
            // case 'H':
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
    if(c=='S' /*|| c=='H' */)
    	endClip = n;
    qStart = beginClip;
    qEnd = beginClip+readLen-1;
    qLen = beginClip+readLen+endClip;
    return contLen;
}

int getAlignmentMatchNum(string cigar, string qSeq, string rSeq, int rStart)
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
    return nMatch;
}

string getCurrentDirectory()
{
	char buff[2000];
	if (getcwd(buff, sizeof(buff)-1) != NULL)
		return buff;
	cerr<< "[getCurrentDirectory] cannot run getcwd()" << endl;
	exit(0);
}

string getExeDirectory()
{
	char buff[2000];
	string str;
	ssize_t len = readlink("/proc/self/exe", buff, sizeof(buff)-1);
	if (len != -1)
	{
		buff[len] = '\0';
		str = buff;
		return str.substr(0, str.find_last_of("/"));
	}
	cerr<< "[getExeDirectory] cannot run readlink()" << endl;
	exit(0);
}