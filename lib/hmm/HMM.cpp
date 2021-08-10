// ***************************************************************************
// HMM.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cmath>
#include <cassert>

#include "HMM.h"
#include "MyDefine.h"
#include "Config.h"
//#include "psiFunc.h"

using namespace std;

void HMM::init() {
	int i, j, k, stateIndx;
	int maxCN = config.getIntPara("maxCN");
	
	// initialize hidden states
	int state_num = 1;
	for(i = 1; i <= maxCN; i++) {
		k = max(i/2, (i+1)/2);
		state_num += i-k+1;
	}
	//cerr << "Number of hidden states in the HMM: " << state_num << endl;
	
	cns.resize(1, state_num, false);
	mafs.resize(1, state_num, false);
	cns.set(0, 0, 0.01);
	mafs.set(0, 0, 0.5);
	stateIndx = 1;
	for(i = 1; i <= maxCN; i++) {
		k = max(i/2, (i+1)/2);
		for(j = k; j <= i; j++) {
			cns.set(0, stateIndx, i);
			mafs.set(0, stateIndx++, 1.0*j/i);
		}
	}
	assert(stateIndx == state_num);
	
	// initialize initial state probability and state transition probability
	pie.resize(1, state_num, false);
	pie.set(1.0/state_num);
	transMat.resize(state_num, state_num, false);
	transMat.set(1);
	double clamp_thres = config.getRealPara("clamp_thres");
	norm_trans(transMat,clamp_thres);
}

double HMM::getObslik_rc(double rc, double lambda, double p) {
	return negabinopdf(rc, lambda, p);
}

double HMM::getObslik_md(double m, double t, double alpha, double beta) {
	return betabinopdf(m, t, alpha, beta);
}

double HMM::ForwardWithScale(Matrix<double> &L, vector<Matrix<double> > &results, int filter) {
	size_t i, j; 	/* state indices */
	size_t t;	/* time index */
	
	int N = L.getROWS(); //numnber of states
	int T = L.getCOLS(); //numnber of observations
	
	results[0].resize(N, T, false); //alpha
	results[2].resize(N, T, false); //gamma
	
	double *pi = pie.getEntrance();
	double *a = transMat.getEntrance();
	double *l = L.getEntrance();
	double *alp = results[0].getEntrance();
	double *gam = results[2].getEntrance();

	double sum;	/* partial sum */
	
	double *scale = new double[T];

	/* 1. Initialization */
	scale[0] = 0.0;
	for(i = 0; i < N; i++) {
		alp[i*T] = pi[i]*l[i*T];
		scale[0] += alp[i*T];
	}
	for(i = 0; i < N; i++) {
		alp[i*T] /= scale[0];
	}
	
	/* 2. Induction */

	for(t = 1; t < T; t++) {
		scale[t] = 0.0;
		for(j = 0; j < N; j++) {
			sum = 0.0;
			for(i = 0; i < N; i++) {
				sum += alp[i*T+t-1]*a[i*N+j];
			}				
			alp[j*T+t]= sum*l[j*T+t];	
			scale[t] += alp[j*T+t];
		}
		for(j = 0; j < N; j++) {
			alp[j*T+t] /= scale[t];
		}
	}

	/* 3. Termination */
	double loglik = 0.0;

	for(t = 0; t < T; t++) {
		loglik += log(scale[t]);
	}
		
	if(filter == 0) {/* forward only */
		for(i = 0; i < N; i++) {
			for(j = 0; j < T; j++) {
				gam[i*T+j] = alp[i*T+j];
			}
		}
	}
	delete[] scale;
	
	return loglik;
}

void HMM::BackwardWithScale(Matrix<double> &L, vector<Matrix<double> > &results) {
	int i, j;   /* state indices */
	int t;      /* time index */
	
	int N = L.getROWS(); //numnber of states
	int T = L.getCOLS(); //numnber of observations
	
	results[1].resize(N, T, false); //beta
	double *bet = results[1].getEntrance();
	
	double *a = transMat.getEntrance();
	double *l = L.getEntrance();
	
	double sum;
	
	/* 1. Initialization */
	
	for(i = 0; i < N; i++) {
		bet[i*T+T-1] = 1.0;
	}
	
	/* 2. Induction */

	for(t = T-2; t >= 0; t--) {
		sum = 0.0;
		for(i = 0; i < N; i++) {
			bet[i*T+t] = 0.0;
			for(j = 0; j < N; j++) {
				bet[i*T+t] += a[i*N+j]*l[j*T+t+1]*bet[j*T+t+1];
			}
			sum += bet[i*T+t];
		}
		for(i = 0; i < N; i++) {
			bet[i*T+t] /= sum;
		}
	}
        	
}

void HMM::ComputeGamma(vector<Matrix<double> > &results) {
/*calculate gamma(:,T) by definition, normalise(alpha(:,T).*beta(:,T))*/
	int i, j;   /* state indices */
	int t;      /* time index */
	
	int N = results[0].getROWS(); //numnber of states
	int T = results[0].getCOLS(); //numnber of observations
	
	//gamma.resize(N,T,false);
	
	double *alp = results[0].getEntrance();
	double *bet = results[1].getEntrance();
	double *gam = results[2].getEntrance();
	
	double sum;

	for(t = 0; t < T; t++) {
		sum = 0.0;
		for(i = 0; i < N; i++) {
			gam[i*T+t] = alp[i*T+t]*bet[i*T+t];
			sum += gam[i*T+t];
		}
		for(i = 0; i < N; i++) {
			gam[i*T+t] /= sum;
		}
	}
		
}

void HMM::ComputeXi(Matrix<double> &L, vector<Matrix<double> > &results) {
	int i, j;
	int t;
	
	int N = results[0].getROWS(); //numnber of states
	int T = results[0].getCOLS(); //numnber of observations
	
	double *a = transMat.getEntrance();
	double *l = L.getEntrance();
	double *alp = results[0].getEntrance();
	double *bet = results[1].getEntrance();
	double *gam = results[2].getEntrance();
	
	Matrix<double> tmp_xi(N, N, false), tmp_xi_summed(N, N, true);
	results[3] = tmp_xi_summed;
	
	double *p1 = tmp_xi.getEntrance();
	double *p2 = results[3].getEntrance();
	
	double sum;
				
	for(t = 0; t <= T-2; t++) {
		sum = 0.0;	
		for(i = 0; i < N; i++) {
			for (j = 0; j < N; j++) {
				p1[i*N+j] = alp[i*T+t]*a[i*N+j]*l[j*T+t+1]*bet[j*T+t+1];
				sum += p1[i*N+j];
			}
		}
		/*calculate gamma with xi for 1 to T-1 (can not calculate gamma(:,T) with this formula) */
		for(i = 0; i < N; i++) {
			gam[i*T+t] = 0.0;
			for(j = 0; j < N; j++) { 
				p1[i*N+j] /= sum;
				p2[i*N+j] += p1[i*N+j];
				gam[i*T+t] += p1[i*N+j];
			}
		}
	}
	
	/*finally calculate gamma(:,T) by definition, normalise(alpha(:,T).*beta(:,T)) */
	sum = 0.0;
	for(i = 0; i < N; i++) {
   		gam[i*T+T-1] = alp[i*T+T-1]*bet[i*T+T-1];
    	sum += gam[i*T+T-1];
	}
	for(i = 0; i < N; i++) {
		gam[i*T+T-1] /= sum;
	}
	
}

double HMM::FwdBack(Matrix<double> &L, vector<Matrix<double> > &results) {
	int i, j;
	int t;
	
	//v[0] = alpha, v[1] = beta, v[2] = gamma, v[3] = xi_summed
	int N = L.getROWS(); //numnber of states
	int T = L.getCOLS(); //numnber of observations
	
	results[0].resize(N, T, false); //alpha
	results[1].resize(N, T, false); //beta
	results[2].resize(N, T, false); //gamma
	
	double *pi = pie.getEntrance();
	double *a = transMat.getEntrance();
	double *l = L.getEntrance();
	double *alp = results[0].getEntrance();
	double *bet = results[1].getEntrance();
	double *gam = results[2].getEntrance();

	double sum, sum1;	/* partial sum */
	
	double *scale = new double[T];

	/* 1. Initialization */
	scale[0] = 0.0;
	for(i = 0; i < N; i++) {
		alp[i*T] = pi[i]*l[i*T];
		scale[0] += alp[i*T];
		bet[i*T+T-1] = 1.0;
	}
	for(i = 0; i < N; i++) {
		alp[i*T] /= scale[0];
	}
	
	/* 2. Induction */

	for(t = 1; t < T; t++) {
		scale[t] = 0.0;
		sum1 = 0.0;
		for(j = 0; j < N; j++) {
			sum = 0.0;
			bet[j*T+T-t-1] = 0.0;
			for(i = 0; i < N; i++) {
				sum += alp[i*T+t-1]*a[i*N+j];
				bet[j*T+T-t-1] += a[j*N+i]*l[i*T+T-t]*bet[i*T+T-t];
			}				
			alp[j*T+t]= sum*l[j*T+t];	
			scale[t] += alp[j*T+t];
			sum1 += bet[j*T+T-t-1];
		}
		for(j = 0; j < N; j++) {
			alp[j*T+t] /= scale[t];
			bet[j*T+T-t-1] /= sum1;
		}
	}

	/* 3. Termination */
	double loglik = 0.0;

	for(t = 0; t < T; t++) {
		loglik += log(scale[t]);
	}
	delete[] scale;
	
	for(t = 0; t < T; t++) {
		sum = 0.0;
		for(i = 0; i < N; i++) {
			gam[i*T+t] = alp[i*T+t]*bet[i*T+t];
			sum += gam[i*T+t];
		}
		for(i = 0; i < N; i++) {
			gam[i*T+t] /= sum;
		}
	}
	
	return loglik;
	
}

double HMM::FwdBack(Matrix<double> &L, vector<Matrix<double> > &results, int filter) {
	//v[0] = alpha, v[1] = beta, v[2] = gamma, v[3] = xi_summed
	int N = L.getROWS(); //numnber of states
	int T = L.getCOLS(); //numnber of observations
	
	/*
	assert(pie.getCOLS() >= pie.getROWS() && pie.getROWS() == 1);
	assert(transMat.getROWS() == transMat.getCOLS());
	assert(pie.getCOLS() == transMat.getROWS());
	assert(N == transMat.getROWS());
	*/
	
	//results[0].resize(N, T, false); //alpha
	//results[2].resize(N, T, false); //gamma
	
	double loglik = ForwardWithScale(L, results, filter);
	
	if(filter > 0) {/*calculate both alpha and beta, but xi_summed is not calculated*/
		BackwardWithScale(L, results);
	}
	
	if(filter == 1) {/*calculate gamma with alpha and beta*/
		ComputeGamma(results);		
	}
	
	if(filter > 1) {/*calculate gamma and xi_summed based on xi*/
		ComputeXi(L, results);
	}
	
	return loglik;
	
}

void HMM::Viterbi(map<string, Matrix<double> > &L_sep, map<string, Matrix<int> > &stateseq_sep) {
	size_t i, j; 	// state indices 
	int t;	// time index 
	int N = pie.getCOLS(); //numnber of states
	
	map<string, Matrix<double> >::iterator it;
	
	double *pi = pie.getEntrance();
	double *a = transMat.getEntrance();
	
	stateseq_sep.clear();
	
	
	for(it = L_sep.begin(); it != L_sep.end(); it++) {
		Matrix<double> &L = (*it).second;
		int T = L.getCOLS(); //numnber of observations
		
		Matrix<double> delta(N, T, false), psi(N, T, false);
		
		double *l = L.getEntrance();
		double *del = delta.getEntrance();
		double *ps = psi.getEntrance();
		
		double sum;
		int maxvalind; 
		double maxval, val;

		/* 1. Initialization */
		
		sum = 0;
		for(i = 0; i < N; i++) {
			del[i*T] = pi[i]*l[i*T];
			ps[i*T] = 0;
			sum += del[i*T];
		}
		for(i = 0; i < N; i++) {
			del[i*T] /= sum;
		}
		
		/* 2. Induction */

		for(t = 1; t < T; t++) {
			sum = 0;
			for(j = 0; j < N; j++) {
				maxval = 0.0;
				maxvalind = 0;
				for(i = 0; i < N; i++) {
					val = del[i*T+t-1]*a[i*N+j];
					if(val > maxval) {
						maxval = val;
						maxvalind = i;
					}
				}				
				del[j*T+t] = maxval*l[j*T+t];
				ps[j*T+t] = maxvalind;
				sum += del[j*T+t];
			}
			for(j = 0; j < N; j++) 
				del[j*T+t] /= sum;
		}

		/* 3. Termination */
		
		Matrix<int> stateseq(1, T, false);
		
		double pprob = 0.0;
		int *q = stateseq.getEntrance();
		q[T-1] = 0;
		
		for(i = 0; i < N; i++) {
			if(del[i*T+T-1] > pprob) {
				pprob = del[i*T+T-1];
				q[T-1] = i;
			}
		}
		
		/* 4. Path (state sequence) backtracking */
		
		for(t = T-2; t >= 0; t--) {
			i = q[t+1];
			q[t] = ps[i*T+t+1];
		}
		
		stateseq_sep.insert(make_pair((*it).first, stateseq));
		
	}
	
}

void HMM::Viterbi(Matrix<double>& L, Matrix<int>& stateseq) {
	size_t i, j; 	/* state indices */
	int t;	/* time index */
	int N = L.getROWS(); //numnber of states
	int T = L.getCOLS(); //numnber of observations
	
	double *pi = pie.getEntrance();
	double *a = transMat.getEntrance();
	
	stateseq.resize(1, T, false);
	
	Matrix<double> delta(N, T, false), psi(N, T, false);
	
	double *l = L.getEntrance();
	double *del = delta.getEntrance();
	double *ps = psi.getEntrance();
	
	double sum;
	int maxvalind; 
	double maxval, val;

	/* 1. Initialization */
	
	sum = 0;
	for(i = 0; i < N; i++) {
		del[i*T] = pi[i]*l[i*T];
		ps[i*T] = 0;
		sum += del[i*T];
	}
	for(i = 0; i < N; i++) {
		del[i*T] /= sum;
	}
	
	/* 2. Induction */

	for(t = 1; t < T; t++) {
		sum = 0;
		for(j = 0; j < N; j++) {
			maxval = 0.0;
			maxvalind = 0;
			for(i = 0; i < N; i++) {
				val = del[i*T+t-1]*a[i*N+j];
				if(val > maxval) {
					maxval = val;
					maxvalind = i;
				}
			}				
			del[j*T+t] = maxval*l[j*T+t];
			ps[j*T+t] = maxvalind;
			sum += del[j*T+t];
		}
		for(j = 0; j < N; j++) 
			del[j*T+t] /= sum;
	}

	/* 3. Termination */
	
	double pprob = 0.0;
	int *q = stateseq.getEntrance();
	q[T-1] = 0;
	
	for(i = 0; i < N; i++) {
		if(del[i*T+T-1] > pprob) {
			pprob = del[i*T+T-1];
			q[T-1] = i;
		}
	}
	
	/* 4. Path (state sequence) backtracking */
	
	for(t = T-2; t >= 0; t--) {
		i = q[t+1];
		q[t] = ps[i*T+t+1];
	}
}

void HMM::assignState(Matrix<double> &gamma, Matrix<int> &stateseq) {
	size_t i; 	/* state indices */
	int t;	/* time index */
	int N = pie.getCOLS(); //number of states
	
	int T = gamma.getCOLS(); //number of observations
	stateseq.resize(1, T, false);
	
	double *p = gamma.getEntrance();
	int *q = stateseq.getEntrance();
	
	for(t = 0; t < T; t++) {
		int maxvalind = 0; 
		double maxval = p[t];
		for(i = 1; i < N; i++) {
			if(p[i*T+t] > maxval) {
				maxval = p[i*T+t];
				maxvalind = i;
			}
		}
		q[t] = maxvalind;
	}
}

void HMM::assignState(map<string, Matrix<double> > &gamma_sep, map<string, Matrix<int> > &stateseq_sep) {
	size_t i; 	/* state indices */
	int t;	/* time index */
	int N = pie.getCOLS(); //numnber of states
	
	map<string, Matrix<double> >::iterator it;
	
	stateseq_sep.clear();
	
	
	for(it = gamma_sep.begin(); it != gamma_sep.end(); it++) {
		Matrix<double> post_probs = (*it).second;
		double *p = post_probs.getEntrance();
		int T = post_probs.getCOLS(); //numnber of observations
		
		Matrix<int> stateseq(1, T, false);
		int *q = stateseq.getEntrance();
		
		for(t = 0; t < T; t++) {
			int maxvalind = 0; 
			double maxval = p[t];
			for(i = 1; i < N; i++) {
				if(p[i*T+t] > maxval) {
					maxval = p[i*T+t];
					maxvalind = i;
				}
			}
			q[t] = maxvalind;
		}
		
		stateseq_sep.insert(make_pair((*it).first, stateseq));		
	}
	
}

