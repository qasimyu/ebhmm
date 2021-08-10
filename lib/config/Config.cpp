// ***************************************************************************
// Config.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

#include "Config.h"

Config::Config() {
	string strParaNames[] = {"ref", "map", "bam", "snp", "input",
							"output", "schmmPath"};
	
	/*---start default configuration---*/
	
	int n = sizeof(strParaNames)/sizeof(string);
	for(int i = 0; i < n; i++) {
		stringParas.insert(make_pair(strParaNames[i], ""));
	}
	
	// parameters for fetching inputs
	intParas.insert(make_pair("baseQ_th", 10));
	intParas.insert(make_pair("mapQ_th", 20));
	intParas.insert(make_pair("minDepth", 5));
	intParas.insert(make_pair("window", 10000));
	intParas.insert(make_pair("male", 0));
	
	// parameters for data preprocessing
	realParas.insert(make_pair("min_td", 5));
	realParas.insert(make_pair("max_td", 300));
	realParas.insert(make_pair("min_gc", 0));
	realParas.insert(make_pair("min_map", 0));
	realParas.insert(make_pair("max_map", 0.98));
	realParas.insert(make_pair("homo_th", 0.01));
	
	// parameters for HMM and EM algorithm
	intParas.insert(make_pair("maxCN", 7));
	intParas.insert(make_pair("ploidy", -1));
	intParas.insert(make_pair("em_max_iter", 30));
	realParas.insert(make_pair("ZERO_FINAL", 2.2204e-16));
	realParas.insert(make_pair("thres_EM", 1e-5));
	realParas.insert(make_pair("thres_del", 0.07));
	realParas.insert(make_pair("clamp_thres", 1-1e-10));
	realParas.insert(make_pair("maxACN", 4.3));
	realParas.insert(make_pair("thres_del", 0.07));
	
	// other parameters
	intParas.insert(make_pair("threads", 1));
	intParas.insert(make_pair("verbose", 1));
	
	/*---end default configuration---*/
}

string Config::getStringPara(string paraName) {
	if(stringParas.find(paraName) != stringParas.end()) {
		return stringParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
	return "";
}

void Config::setStringPara(string paraName, string value) {
	if(stringParas.find(paraName) != stringParas.end()) {
		stringParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

long Config::getIntPara(string paraName) {
	if(intParas.find(paraName) != intParas.end()) {
		return intParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setIntPara(string paraName, long value) {
	if(intParas.find(paraName) != intParas.end()) {
		intParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}

double Config::getRealPara(string paraName) {
	if(realParas.find(paraName) != realParas.end()) {
		return realParas[paraName];
	}
	else {
		cerr << "Error: unrecognized parameter name \"" << paraName << "\"" << endl;
		exit(1);
	}
}

void Config::setRealPara(string paraName, double value) {
	if(realParas.find(paraName) != realParas.end()) {
		realParas[paraName] = value;
	}
	else {
		cerr << "Warning: unrecognized parameter name \"" << paraName << "\"" << endl;
	}
}


