// ***************************************************************************
// GenomeData.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _GENOMEDATA_H
#define _GENOMEDATA_H

#include <vector>
#include <map>
#include <cstring>
#include <pthread.h>

#include "Matrix.h"

class GenomeData {
	private:
		string sampleName; //the name of sample to be analyzed
		vector<int> chromosomes; //chromosomes names
		/*
		map<int, Matrix<double> > data_pos_sep; //snp positions
		map<int, Matrix<double> > data_rc_sep; //read counts
		map<int, Matrix<double> > data_md_sep; //major allele read depth
		map<int, Matrix<double> > data_td_sep; //total allele read depth
		*/
		map<int, Matrix<double> > data_sep; // data grouped by chr: positions, major allele and total read depths, read counts
		map<int, Matrix<double> > data_sep_ds;
		map<int, Matrix<double> > state_segments; // state segments
		
		int num_snp;
		double rc_median, max_rc;
		int stepsize_ds;
		
		pthread_mutex_t pm;
		
		void corrRC(double *rcs_p, double *factors, int T);
		void preprocessRD(double *bd_p, double *td_p, int T);
		void filterSNPs(Matrix<double> &data_all);
		void splitDataByChr(Matrix<double> &data_all);
		void saveData(Matrix<double> &data_all);
		
		int readSNPs(map<string, vector<long long> > &positions_all, vector<string> &chroms);
		void calculateGC(string sequence, double &gc_content, long &length);
		void calculateMapScore(double *chromData, long long spos, long long epos, double &mapScore);
		
		static void* smoothRCByChr(const void *arg);
		
	public:
		GenomeData() {pthread_mutex_init(&pm, NULL);}
		~GenomeData() {}
		
		void LoadAndProcess();
		void saveResults();
		
		string getSampleName() {return sampleName;}
		vector<int>& getChroms() {return chromosomes;}
		int getChrCount() {return chromosomes.size();}
		
		double getMedRC() {return rc_median;}
		double getMaxRC() {return max_rc;}
		
		int getStepSizeDS() {return stepsize_ds;}
		void setStepSizeDS(int n) {stepsize_ds = n;}
		
		double calculateACN();
		
		Matrix<double>& getChrData(int chr);
		
		Matrix<double>& getSegments(int chr) {return state_segments[chr];}
		void setSegments(int chr, Matrix<double> &segments);
		
		void fetchInputs();
};


#endif

