// ***************************************************************************
// InputParser.cpp (c) 2019 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <string>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#include <limits.h>

#include "InputParser.h"
#include "MyDefine.h"

using namespace std;

void InputParser::parseArgs(int argc, char *argv[]) {
	if(argc == 1) {
		usage(argv[0]);
		exit(0);
	}
	string subcmd = argv[1];
	if(subcmd.compare("-h") == 0 || subcmd.compare("--help") == 0) {
		usage(argv[0]);
		exit(0);
	}
	else if(subcmd.compare("-v") == 0 || subcmd.compare("--version") == 0) {
		cerr << "EBHMM version " << current_version << endl;
		exit(0);
	}
	else if(subcmd.compare("prepareinputs") == 0) {
		parseArgs_prepareInputs(argc-1, &argv[1]);
	}
	else if(subcmd.compare("callcnas") == 0) {
		parseArgs_callCNAs(argc-1, &argv[1]);
	}
	else {
		cerr << "Error: unrecognized subcommand \"" << subcmd << "\"." << endl;
		usage(argv[0]);
		exit(0);
	}
	
	string binaryPath = argv[0];
	char abs_path_buff[PATH_MAX];
	realpath(binaryPath.c_str(), abs_path_buff);
	string ebhmmPath = abs_path_buff;
	int indx = ebhmmPath.rfind("bin");
	if(indx != string::npos) {
		ebhmmPath = ebhmmPath.substr(0, indx-1);
	}
	else {
		indx = ebhmmPath.rfind("/");
		ebhmmPath = ebhmmPath.substr(0, indx);
	}
	//cerr << ebhmmPath << endl;
	config.setStringPara("ebhmmPath", ebhmmPath);
}

void InputParser::parseArgs_prepareInputs(int argc, char *argv[]) {
	string refFile = "", bwFile = "";
	string bamFile = "", snpFile = "", outFile = "";
	int window = 10000, ismale = 0;
	int baseQ_th = 10, mapQ_th = 20, minDepth = 5;

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"bam", required_argument, 0, 'b'},
		{"ref", required_argument, 0, 'r'},
		{"map", required_argument, 0, 'm'},
		{"snp", required_argument, 0, 's'},
		{"male", no_argument, 0, 'M'},
		{"window", required_argument, 0, 'w'},
		{"baseQ", required_argument, 0, 'Q'},
		{"mapQ", required_argument, 0, 'q'},
		{"minDepth", required_argument, 0, 'd'},
		{"output", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hb:r:m:s:Mw:Q:q:d:o:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage_prepareInputs(argv[0]);
				exit(0);
			case 'b':
				bamFile = optarg;
				break;
			case 'r':
				refFile = optarg;
				break;
			case 'm':
				bwFile = optarg;
				break;
			case 's':
				snpFile = optarg;
				break;
			case 'M':
				ismale = 1;
				break;
			case 'w':
				window = atoi(optarg);
				break;
			case 'Q':
				baseQ_th = atoi(optarg);
				break;
			case 'q':
				mapQ_th = atoi(optarg);
				break;
			case 'd':
				minDepth = atoi(optarg);
				break;
			case 'o':
				outFile = optarg;
				break;
			default :
				usage_prepareInputs(argv[0]);
				exit(1);
		}
	}
	
	if(bamFile.empty()){
        cerr << "Use --bam to specify BAM file (.bam)." << endl;
		usage_prepareInputs(argv[0]);
        exit(1);
    }

	if(refFile.empty()){
		cerr << "Use --ref to specify reference file (.fasta)." << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	if(bwFile.empty()){
		cerr << "Use --map to specify mappability file(.bigwig)." << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	if(snpFile.empty()){
		cerr << "Use --snp to specify SNP file." << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	if(outFile.empty()){
		cerr << "Use --output to specify the output file." << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	if(window < 1000 || window > 500000) {
		cerr << "Error: the value of parameter \"window\" should be in 1000~500000!" << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	if(ismale == 0) {
		cerr << "Chromosome Y will not be analyzed." << endl;
	}
	
	if(baseQ_th < 0) {
		cerr << "Error: the value of parameter \"baseQ\" should be a nonnegative integer!" << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	if(mapQ_th < 0) {
		cerr << "Error: the value of parameter \"mapQ\" should be a nonnegative integer!" << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	if(minDepth < 0) {
		cerr << "Error: the value of parameter \"minDepth\" should be a nonnegative integer!" << endl;
		usage_prepareInputs(argv[0]);
		exit(1);
	}
	
	config.setStringPara("ref", refFile);
	config.setStringPara("map", bwFile);
	config.setStringPara("bam", bamFile);
	config.setStringPara("snp", snpFile);
	config.setStringPara("output", outFile);
	config.setIntPara("window", window);
	config.setIntPara("male", ismale);
	config.setIntPara("mapQ_th", mapQ_th);
	config.setIntPara("baseQ_th", baseQ_th);
	config.setIntPara("minDepth", minDepth);
}

void InputParser::parseArgs_callCNAs(int argc, char *argv[]) {
	string dataFile = "", outputDir = "";
	int minDepth = 5, maxDepth = 300;
	int maxCN = 7, ploidy = -1;
	int threads = 1;

	struct option long_options[] = {
		{"help", no_argument, 0, 'h'},
		{"input", required_argument, 0, 'i'},
		{"minDepth", required_argument, 0, 'd'},
		{"maxDepth", required_argument, 0, 'D'},
		{"maxCN", required_argument, 0, 'C'},
		{"ploidy", required_argument, 0, 'p'},
		{"threads", required_argument, 0, 't'},
		{"output", required_argument, 0, 'o'},
		{0, 0, 0, 0}
	};

	int c;
	//Parse command line parameters
	while((c = getopt_long(argc, argv, "hi:d:D:C:p:t:o:", long_options, NULL)) != -1){
		switch(c){
			case 'h':
				usage_callCNAs(argv[0]);
				exit(0);
			case 'i':
				dataFile = optarg;
				break;
			case 'd':
				minDepth = atoi(optarg);
				break;
			case 'D':
				maxDepth = atoi(optarg);
				break;
			case 'C':
				maxCN = atoi(optarg);
				break;
			case 'p':
				ploidy = atoi(optarg);
				break;
			case 't':
				threads = atoi(optarg);
				break;
			case 'o':
				outputDir = optarg;
				break;
			default :
				usage_callCNAs(argv[0]);
				exit(1);
		}
	}

	if(dataFile.empty()){
        cerr << "Use --input to specify input file." << endl;
		usage_callCNAs(argv[0]);
        exit(1);
    }
	
	if(minDepth < 0) {
		cerr << "Error: the value of parameter \"minDepth\" should be a nonnegative integer!" << endl;
		usage_callCNAs(argv[0]);
		exit(1);
	}
	
	if(maxDepth < minDepth) {
		cerr << "Error: the value of parameter \"maxDepth\" should be larger than \"minDepth\"!" << endl;
		usage_callCNAs(argv[0]);
		exit(1);
	}
	
	if(maxCN < 3) {
		cerr << "Error: the value of parameter \"maxCN\" should be larger than 2!" << endl;
		usage_callCNAs(argv[0]);
		exit(1);
	}
	
	if(ploidy == -1) {
		cerr << "ploidy will be inferred from the data." << endl;
	}
	else if(ploidy != 2 && ploidy != 3 && ploidy != 4) {
		cerr << "Error: the value of parameter \"ploidy\" should be in [2,3,4]!" << endl;
		usage_callCNAs(argv[0]);
		exit(1);
	}
	else {
		cerr << "ploidy is set to " << ploidy << endl;
	}
	
	if(threads < 1) {
		cerr << "Error: the value of parameter \"threads\" should be a positive integer!" << endl;
		usage_callCNAs(argv[0]);
		exit(1);
	}
	
	if(outputDir.empty()){
		cerr << "Use --output to specify the output directory." << endl;
		usage_callCNAs(argv[0]);
		exit(1);
	}
	
	bool is_exist = (access(outputDir.c_str(), F_OK) == 0);
	if(!is_exist) {
		cerr << "Error: the output directory " << outputDir << " does not exist!" << endl;
		exit(1);
	}
	//string cmd = "test ! -e "+outputDir+" && mkdir -m 755 -p "+outputDir;
	//system(cmd.c_str());
	
	mkdir((outputDir+"/plots").c_str(), 0755);
	
	/*** create thread pool ***/
	threadpool = new ThreadPool(threads);
	threadpool->pool_init();
	
	config.setStringPara("input", dataFile);
	config.setStringPara("output", outputDir);
	config.setRealPara("min_td", minDepth);
	config.setRealPara("max_td", maxDepth);
	config.setIntPara("maxCN", maxCN);
	config.setIntPara("ploidy", ploidy);
	config.setIntPara("threads", threads);
}

void InputParser::usage(const char* app) {
	cerr << "\nEBHMM version: " << current_version << endl;
	cerr << "Usage: " << app << " [subcommand] [options]" << endl
		<< endl
		<< "Optional arguments:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -v, --version <string>          print software version" << endl
		<< endl
		<< "Available subcmds:" << endl
		<< "    prepareinputs        fetch inputs from tumor sequencing data for calling CNAs" << endl
		<< "    callcnas             detect copy number alterations" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

void InputParser::usage_prepareInputs(const char* app) {
	cerr << "Usage: ebhmm " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -b, --bam <string>              BAM file" << endl
		<< "    -r, --ref <string>              genome reference file(.fasta)" << endl
		<< "    -m, --map <string>              mappability file(.bw)" << endl
		<< "    -s, --snp <string>              SNP file" << endl
		<< "    -M, --male                      indicate sample from a male" << endl
		<< "    -w, --window <int>              set the size of windows [default:10000]" << endl
		<< "    -Q, --baseQ <int>               threshold value for base quality [default:10]" << endl
		<< "    -q, --mapQ <int>                threshold value for mapping quality [default:20]" << endl
		<< "    -d, --minDepth <int>            minimum read depth for a position to be considered [default:5]" << endl
		<< "    -o, --output <string>           output file" << endl
		<< endl
		<< "Example:" << endl
		<< "    ebhmm " << app << " -b example.bam -r hg19.fasta -m hg19.bw -s hg19.snp151.txt -M -w 50000 -o example.out" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}

void InputParser::usage_callCNAs(const char* app) {
	cerr << "Usage: ebhmm " << app << " [options]" << endl
		<< endl
		<< "Options:" << endl
		<< "    -h, --help                      give this information" << endl
		<< "    -i, --input <string>            input file generated by \"prepareInputs\" program " << endl
		<< "    -d, --minDepth <int>            filter out positions with read depth less than this value [default:5]" << endl
		<< "    -D, --maxDepth <int>            filter out positions with read depth larger than this value [default:300]" << endl
		<< "    -C, --maxCN <int>               maximum copy number to be considered [default:7]" << endl
		<< "    -p, --ploidy <int>              ploidy of the sample, should be a value in (2, 3, 4) [default:to be inferred from the data]" << endl
		<< "    -t, --threads <int>             number of threads to use [default:1]" << endl
		<< "    -o, --output <string>           output directory to save results" << endl
		<< endl
		<< "Example:" << endl
		<< "    ebhmm " << app << " -i example.out -d 10 -D 200 -C 5 -p 2 -t 5 -o results" << endl
		<< endl
		<< "Author: Zhenhua Yu <qasim0208@163.com>\n" << endl;
}
