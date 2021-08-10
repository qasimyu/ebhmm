// ***************************************************************************
// Model.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _MODEL_H
#define _MODEL_H

#include <vector>
#include <map>
#include <pthread.h>

#include "HMM.h"
#include "Matrix.h"

class ModelParas {
	public:
		Matrix<double> mu; //expected major allele ratio
		Matrix<double> eta; //measure of variance in mu
		double lambda; //mean read counts of copy neutral regions
		Matrix<double> p; //parameters of NBD
		Matrix<double> normal_prior; //parameters for observation function, prior(0,2,4,6 copies)
		Matrix<int> indicator; //indicator vector: '1' for update, '0' for fixed
};

class Model {
	private:
		int modelIndx;
		HMM hmm;
		ModelParas paras, n_paras;	
		Matrix<double> lambda_limit;
		
		pthread_mutex_t pm;
		
		static int nr_max_iter;
		
		double loglik;
		map<int, Matrix<double> > gamma_sep;
		map<int, Matrix<double> > condi_probs_sep;
		map<int, Matrix<double> > condi_probs_fluct_sep;
		Matrix<double> exp_num_visits;
		Matrix<double> exp_num_trans;
		
		static void* fwdBackPerChr(const void *arg);
		static void* update_mu(const void *arg);
		static void* update_lambda(const void *arg);
		static void* update_eta(const void *arg);
		static void* update_p(const void *arg);
		static void* predictChrStates(const void *arg);
		
		void getObslik(Matrix<double> &rd, Matrix<double> &rc, vector<Matrix<double> > &results);
		void getObslik(int T, double *rd, double *rc, vector<Matrix<double> > &results);
		void getSegments(Matrix<int> &stateseq, Matrix<double> &segments);
		void correctStates(int chr, Matrix<double> &segments);
		void calculateScores(int chr, Matrix<double> &segments);

	public:
		void initParas(int modelIndx, double lambda);
		double estimateParas();
		void printParas();
		void predictStates();
		static bool em_converged(double ll, double pre_ll);
		
		void compute_ess();
		
		void process_results(Matrix<double> &p_states, double &acn, map<int, Matrix<double> > &segments_sep);
		
		HMM& getHMM() {return hmm;}
		ModelParas& getParas() {return paras;}

		pthread_mutex_t& getPM() {return pm;}
		
};


#endif

