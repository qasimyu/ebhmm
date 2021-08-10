// ***************************************************************************
// CNACaller.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "CNACaller.h"
#include "MyDefine.h"
#include "Config.h"
#include "Matrix.h"
//#include "MathFunc.h"

using namespace std;

void CNACaller::callCNAs() {
	init();
	getBestModel();
	predict();
	plotResults();
	//plotResults_new();
}

void CNACaller::init() {
	int i;
	double rc_med = genomedata.getMedRC();
	double lambda;
	
	vector<int> ploidy_all;
	int ploidy = config.getIntPara("ploidy");
	if(ploidy != -1) {
		ploidy_all.push_back(ploidy);
	}
	else {
		ploidy_all.push_back(2); //diploidy
		ploidy_all.push_back(3); //triploidy
		ploidy_all.push_back(4); //tetraploidy
	}
	num_model = ploidy_all.size();
	models = new Model[num_model];
	Model baseModel;
	for(i = 0; i < num_model; i++) {
		ploidy = ploidy_all[i];
		lambda = 2*rc_med/ploidy;
		models[i] = baseModel;
		models[i].initParas(i, lambda);
	}
	LL_all.resize(1, num_model, false);
	acn_all.resize(1, num_model, false);
	
	int stateCount = models[0].getHMM().getStateCount();
	p_states_all.resize(stateCount, num_model, false);
	
	p_total_del.resize(1, num_model, false);
	valid_indic.resize(1, num_model, false);
	valid_indic.set(1);
	
	pthread_mutex_init(&pm, NULL);
}

void CNACaller::getBestModel() {
	int i, j, k;
	
	//search for lambda values
	cerr << "\nfinding optimal model parameters...\n" << endl;
	for(i = 0; i < num_model; i++) {
		inferModelParas(i);
	}
	
	//correct for bias that higher copy numbers are more preferred due to more candidate genotypes
	int num_snp = 0;
	vector<int>& chroms = genomedata.getChroms();
	for(i = 0; i < chroms.size(); i++) {
		num_snp += genomedata.getChrData(chroms[i]).getCOLS();
	}
	Matrix<double> &cns = models[0].getHMM().getCNs();
	Matrix<double> tmp(1, cns.getCOLS());
	for(i = 0; i < cns.getCOLS(); i++) {
		k = ceil((cns[i]+0.01)/2);
		k = cns[i]-k+2;
		tmp[i] = 1-1.0/k;
	}
	Matrix<double> df_bias = tmp*p_states_all;
	
	//check solution validity
	noSolutionFlag = false;
	if(valid_indic.sum() == 0) {
		noSolutionFlag = true;
		valid_indic.set(1);
		cerr << "Warning: cannot find a feasible solution with pre-defined criteria!" << endl;
	}
	
	//fetch the solution that have highest LL
	double best_ll = -1.0e15;
	for(i = 0; i < num_model; i++) {
		if(valid_indic[i] == 0) {
			continue;
		}
		if(LL_all[i] > best_ll) {
			best_ll = LL_all[i];
			best_indx = i;
		}
	}
	
	//correct for df bias
	double thres_ratio = 0.0008;
	Matrix<double> ratios(1, num_model);
	double candi_best = best_indx;
	double best_ratio = -1.0e10;
	for(i = 0; i < num_model; i++) {
		if(valid_indic[i] == 0) {
			continue;
		}
		double tmp1 = best_ll-LL_all[i];
		double tmp2 = log(num_snp)/2*(df_bias[candi_best]-df_bias[i]);
		if(fabs(tmp1) <= 1 && fabs(tmp2) <= 0.01) {
			ratios[i] = 0.0;
		}
		else {
			ratios[i] = tmp2/tmp1;
		}
		
		if(ratios[i] >= thres_ratio && ratios[i] > best_ratio) {
			best_ratio = ratios[i];
			best_indx = i;
		}
	}
	
	/*
	// for debugging, save intermediate results
	ofstream ofs;
	string outputFile = "model.details";
	ofs.open(outputFile.c_str(), ofstream::out | ofstream::app);
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << outputFile << endl;
		exit(-1);
	}
	ofs << "----------------------------------------------------------" << endl;
	ofs << genomedata.getSampleName() << endl;
	for(i = 0; i < num_model; i++) {
		if(i == best_indx) {
			ofs << 2;
		}
		else if(i == candi_best) {
			ofs << 1;
		}
		else if(valid_indic.get(0,i) == 1) {
			ofs << 0;
		}
		else {
			ofs << -1;
		}
		ofs << "\t" << p_total_del[i] << "\t" << acn_all[i] << "\t" << ratios[i] << "\t" << LL_all[i] << endl;
	}
	ofs << "----------------------------------------------------------" << endl;
	ofs.close();
	*/
}

void CNACaller::inferModelParas(int mindx) {
	double thres_del = config.getRealPara("thres_del");
	double maxACN = config.getRealPara("maxACN");
	
	Matrix<double> p_states;
	double acn;
	map<int, Matrix<double> > segments_sep;
	
	double loglik = models[mindx].estimateParas();
	LL_all[mindx] = loglik;
	models[mindx].process_results(p_states, acn, segments_sep);
	p_states_all.setCol(mindx, p_states);
	acn_all[mindx] = acn;
	p_total_del[mindx] = p_states.get(0,0);
	if(p_states.get(0,0) > thres_del || acn > maxACN) {
		valid_indic[mindx] = 0;
	}
	
	cerr << endl;
	cerr << "--------------- screening report -----------------" << endl;
	cerr << "run " << mindx+1 << " done." << endl;
	models[mindx].printParas();
	//cerr << "LL: " << loglik << endl;
	printf("LL: %.2f\n", loglik);
	cerr << "--------------- screening report -----------------" << endl;
}

void CNACaller::predict() {
	int i, j;
	
	//use best model to fine-tune paras and call CNAs
	genomedata.setStepSizeDS(1);
	curModel_p = &models[best_indx];
	models[best_indx].estimateParas();
	cerr << "\nInference of optimal model parameters done!" << endl;
	
	models[best_indx].predictStates();
	double acn = genomedata.calculateACN();
	
	string sampleName = genomedata.getSampleName();
	
	//write brief information to log file
	ofstream ofs;
	string outputFile = "log.txt";
	struct stat buffer;
	bool is_exist = (stat(outputFile.c_str(), &buffer) == 0);
	ofs.open(outputFile.c_str(), ofstream::out | ofstream::app);
	if(!ofs.is_open()) {
		cerr << "Warning: cannot open file " << outputFile << endl;
	}
	else {
		if(!is_exist) {
			ofs << "#Version\tDate\tTime\tSample\tAverage_copy_number" << endl;
		}
		ofs << "SCHMM-" << current_version;
		struct timeval tv;
		char buf[64];
		gettimeofday(&tv, NULL);
		strftime(buf, sizeof(buf)-1, "%Y-%m-%d\t%H:%M:%S", localtime(&tv.tv_sec));
		ofs << "\t" << buf << "\t" << sampleName << "\t" << acn << endl;
	}
	ofs.close();
	
	ModelParas &paras = models[best_indx].getParas();
	//save CNA detection results
	string outputDir = config.getStringPara("output");
	outputFile = outputDir+"/"+sampleName+".result";
	FILE *fp = fopen(outputFile.c_str(),"w");
	if(!fp) {
		cerr << "Error: cannot open file " << outputFile << " to save the results!" << endl;
		exit(-1);
	}
	fprintf(fp, "---------------------------------------------------------------\n");
	fprintf(fp, "\t\t\tSummary of SCHMM results (version %s)\n", current_version.c_str());
	if(noSolutionFlag) {
		fprintf(fp, "Warning: Prediction results may be inaccurate due to the failure\n");
		fprintf(fp, "in finding feasible model parameters!\n");
	}
	fprintf(fp, "General information of this tumor sample:\n");
	fprintf(fp, "\tAverage copy number: %.2f\n", acn);
	fprintf(fp, "\tCopy neutral read counts: %.2f\n", paras.lambda);
	fprintf(fp, "\tMus of BBDs for major allele depth:\n\t");
	fprintf(fp, "%.6f", paras.mu.get(0, 0));
	for(i = 1; i < paras.mu.getCOLS(); i++) {
		fprintf(fp, " %.6f", paras.mu.get(0, i));
	}
	fprintf(fp, "\n");
	fprintf(fp, "\tEtas of BBDs for major allele depth:\n\t");
	fprintf(fp, "%.6f", paras.eta.get(0, 0));
	for(i = 1; i < paras.eta.getCOLS(); i++) {
		fprintf(fp, " %.6f", paras.eta.get(0, i));
	}
	fprintf(fp, "\n");
	fprintf(fp, "\tParameter p of NB distributions:\n\t");
	fprintf(fp, "%.6f", paras.p.get(0, 0));
	for(i = 1; i < paras.p.getCOLS(); i++) {
		fprintf(fp, " %.6f", paras.p.get(0, i));
	}
	fprintf(fp, "\n");
	fprintf(fp, "---------------------------------------------------------------\n");
	fprintf(fp, "##Aberration=<Value=HOMD,Description=\"homozygous deletion\">\n");
	fprintf(fp, "##Aberration=<Value=HEMD,Description=\"hemizygous deletion\">\n");
	fprintf(fp, "##Aberration=<Value=NHET,Description=\"copy neutral heterozygosity\">\n");
	fprintf(fp, "##Aberration=<Value=NLOH,Description=\"copy neutral LOH\">\n");
	fprintf(fp, "##Aberration=<Value=ALOH,Description=\"copy amplified LOH\">\n");
	fprintf(fp, "##Aberration=<Value=AHET,Description=\"copy amplified heterozygosity\">\n");
	fprintf(fp, "#Chr\tStartPos\tEndPos\tCN\tmCN\tAberration\tScore\n");
	vector<int>& chroms = genomedata.getChroms();
	int chr, num_seg;
	int sindx, eindx, cn, mcn;
	double score;
	string aber_type;
	for(i = 0; i < chroms.size(); i++) {
		chr = chroms[i];
		Matrix<double> &chrData = genomedata.getChrData(chr);
		double *p1 = chrData.getEntrance();
		
		Matrix<double> &segs = genomedata.getSegments(chr);
		num_seg = segs.getROWS();
		for(j = 0; j < num_seg; j++) {
			sindx = segs.get(j, 1);
			eindx = segs.get(j, 2);
			cn = segs.get(j, 3);
			mcn = segs.get(j, 4);
			score = segs.get(j, 5);
			
			if(cn == 0) {
				aber_type = "HOMD";
			}
			else if(cn == 1) {
				aber_type = "HEMD";
			}
			else if(cn == 2 && mcn == 1) {
				aber_type = "NHET";
			}
			else if(cn == 2 && mcn == 2) {
				aber_type = "NLOH";
			}
			else if(cn == mcn) {
				aber_type = "ALOH";
			}
			else {
				aber_type = "AHET";
			}
			
			fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%s\t%.2f\n", chr, (int) p1[sindx], (int) p1[eindx], cn, mcn, aber_type.c_str(), score);
		}
	}
	fclose(fp);
}

/*
void CNACaller::plotResults() {
	string sampleName = genomedata.getSampleName();
	string outputDir = config.getStringPara("output");
	string scriptFile = outputDir+"/"+sampleName+".m";
	
	ofstream ofs;
	ofs.open(scriptFile.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: can not open file " << scriptFile << endl;
		exit(-1);
	}
    string dataFile = outputDir+"/"+sampleName+".processed.txt";
    string resultFile = outputDir+"/"+sampleName+".result";
    string plotDir = outputDir+"/plots";
	
	string schmmPath = config.getStringPara("schmmPath");
	int maxCN = config.getIntPara("maxCN");
	ofs << "path(path, '" << schmmPath << "/util/matlab_code')" << endl;
    ofs << "plotResults(" << maxCN << ", '" << dataFile << "', '" << resultFile << "', '" << plotDir << "', '" << sampleName << "')";
	
	ofs.close();
	
	cerr << "\nplot results now..." << endl;
	
	//string cmd = "matlab -nodisplay -nosplash < "+scriptFile+" 1>/dev/null 2>"+outputDir+"/running.err";
	string cmd = "matlab -nodisplay -nosplash < "+scriptFile+" 1>/dev/null 2>/dev/null";
	system(cmd.c_str());
	
	remove(scriptFile.c_str());
}
*/

void CNACaller::plotResults() {
	cerr << "\nplot results now..." << endl;
	vector<int>& chroms = genomedata.getChroms();
	for(int i = 0; i < chroms.size(); i++) {
		plotChrResults(chroms[i]);
	}
	plotChrResults(0);
}

void CNACaller::plotResults_new() {
	cerr << "\nplot results now..." << endl;
	vector<int>& chroms = genomedata.getChroms();
	for(int i = 0; i < chroms.size(); i++) {
		threadpool->pool_add_work(&CNACaller::plotChrResults, &chroms[i], i);
	}
	threadpool->pool_add_work(&CNACaller::plotChrResults, NULL, chroms.size());
	threadpool->wait();
}

void CNACaller::plotChrResults(int chr) {
	string sampleName = genomedata.getSampleName();
	string outputDir = config.getStringPara("output");
	char scriptFile[200];
	sprintf(scriptFile, "%s/%s_%d.m", outputDir.c_str(), sampleName.c_str(), chr);
	
	string dataFile = outputDir+"/"+sampleName+".processed.txt";
    string resultFile = outputDir+"/"+sampleName+".result";
    string plotDir = outputDir+"/plots";
	string schmmPath = config.getStringPara("schmmPath");
	int maxCN = config.getIntPara("maxCN");
	
	ofstream ofs;
	ofs.open(scriptFile);
	if(!ofs.is_open()) {
		cerr << "Error: can not open file " << scriptFile << endl;
		exit(-1);
	}
	ofs << "path(path, '" << schmmPath << "/src')" << endl;
	//ofs << "cd " << schmmPath << "/src" << endl;
    ofs << "plotResults(" << chr << ", " << maxCN << ", '" << dataFile << "', '" << resultFile << "', '" << plotDir << "', '" << sampleName << "')";
	ofs.close();
	
	string cmd = "matlab -nodisplay -nosplash < "+string(scriptFile)+" 1>/dev/null 2>/dev/null";
	system(cmd.c_str());
	remove(scriptFile);
}

void* CNACaller::plotChrResults(const void *arg) {
	int chr;
	if(arg == NULL) {
		chr = 0;
	}
	else {
		chr = *((int*) arg);
	}
	
	string sampleName = genomedata.getSampleName();
	string outputDir = config.getStringPara("output");
	char scriptFile[200];
	sprintf(scriptFile, "%s/%s_%d.m", outputDir.c_str(), sampleName.c_str(), chr);
	
	string dataFile = outputDir+"/"+sampleName+".processed.txt";
    string resultFile = outputDir+"/"+sampleName+".result";
    string plotDir = outputDir+"/plots";
	string schmmPath = config.getStringPara("schmmPath");
	int maxCN = config.getIntPara("maxCN");
	
	ofstream ofs;
	ofs.open(scriptFile);
	if(!ofs.is_open()) {
		cerr << "Error: can not open file " << scriptFile << endl;
		exit(-1);
	}
	ofs << "path(path, '" << schmmPath << "/src')" << endl;
	//ofs << "cd " << schmmPath << "/src" << endl;
    ofs << "plotResults(" << chr << ", " << maxCN << ", '" << dataFile << "', '" << resultFile << "', '" << plotDir << "', '" << sampleName << "')";
	ofs.close();
	
	string cmd = "matlab -nodisplay -nosplash < "+string(scriptFile)+" 1>/dev/null 2>/dev/null";
	system(cmd.c_str());
	remove(scriptFile);
}

