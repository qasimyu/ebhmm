// ***************************************************************************
// HMM.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _HMM_H
#define _HMM_H

#include <map>
#include "Matrix.h"

class HMM {
	private:
		//hidden states
		Matrix<double> cns;
		Matrix<double> mafs;
		
		//pie and trans
		Matrix<double> pie;
		Matrix<double> transMat;
		
		/* functions that will be used by FwdBack algorithm */
		double ForwardWithScale(Matrix<double> &L, vector<Matrix<double> > &results, int filter);
		void BackwardWithScale(Matrix<double> &L, vector<Matrix<double> > &results);
		void ComputeGamma(vector<Matrix<double> > &results);			
		void ComputeXi(Matrix<double> &L, vector<Matrix<double> > &results);
		
	public:
		HMM() {}
		void init();
		void resetProbs();
		
		double getObslik_rc(double rc, double lambda, double p); // read counts, lambda, p
		double getObslik_md(double m, double t, double alpha, double beta); // major-allele depth, total depth, alpha, beta
		
		double FwdBack(Matrix<double> &L, vector<Matrix<double> > &results);
		double FwdBack(Matrix<double> &L, vector<Matrix<double> > &results, int filter);
		
		void Viterbi(Matrix<double> &L, Matrix<int> &stateseq);
		void Viterbi(map<string, Matrix<double> > &L_sep, map<string, Matrix<int> > &stateseq_sep);
		void assignState(map<string, Matrix<double> > &gamma_sep, map<string, Matrix<int> > &stateseq_sep);
		void assignState(Matrix<double> &gamma, Matrix<int> &stateseq);
		
		Matrix<double>& getCNs() {return cns;}
		Matrix<double>& getMAFs() {return mafs;}
		
		Matrix<double>& getPI() {return pie;}
		Matrix<double>& getTransMat() {return transMat;}
		void setPI(Matrix<double> &pie) {this->pie = pie;}
		void setTransMat(Matrix<double> &transMat) {this->transMat = transMat;}
		int getStateCount() {return pie.getCOLS();}
};


#endif

