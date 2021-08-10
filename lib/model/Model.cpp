// ***************************************************************************
// Model.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <cmath>

#include "Model.h"
#include "MyDefine.h"

using namespace std;

int Model::nr_max_iter = 10;

void Model::initParas(int modelIndx, double lambda) {
	pthread_mutex_init(&pm, NULL);
	this->modelIndx = modelIndx;
	
	int i;
	
	//initialize HMM
	hmm.init();
	
	Matrix<double>& mafs = hmm.getMAFs();
	int num_state = mafs.getCOLS();
	
	exp_num_visits.resize(1, num_state, false);
	exp_num_trans.resize(num_state, num_state, false);
	
	//initialize mu
	paras.mu.resize(1, num_state, false);
	for(i = 0; i < num_state; i++) {
		if(fabs(mafs.get(0, i)-1) < 1e-5) {
			paras.mu.set(0, i, 0.97);
		}
		else {
			paras.mu.set(0, i, mafs.get(0, i));
		}
	}
	
	//initialize eta
	paras.eta.resize(1, num_state, false);
	paras.eta.set(0.01);
	
	//initialize lambda
	paras.lambda = lambda;
	lambda_limit.resize(1, 2, false);
	lambda_limit.set(0, 0, 0.8*paras.lambda);
	lambda_limit.set(0, 1, 1.2*paras.lambda);
	
	//initialize p
	int maxCN = config.getIntPara("maxCN");
	paras.p.resize(1, maxCN+1, false);
	paras.p.set(0.9);
	
	//initialize normal_prior
	paras.normal_prior.resize(1, 1+maxCN/2, false);
	paras.normal_prior.set(0, 0, 1);
	double med_rc = genomedata.getMedRC();
	for(i = 1; i < 1+maxCN/2; i++) {
		double v;
		if(fabs(lambda/med_rc-1) < 0.01) { //diploidy
			v = 1.7;
		}
		else if(fabs(lambda/med_rc-0.5) < 0.01) { //tetraploidy
			v = 1.65;
		}
		else { // triploidy
			v = 1.7;
		}
		
		if(v < 1.0) {
			paras.normal_prior.set(0, i, 1.0);
		}
		else {
			paras.normal_prior.set(0, i, v);
		}
	}
	
	//initialize indicator
	paras.indicator.resize(1, 6, false);
	paras.indicator.set(1);
}

bool Model::em_converged(double ll, double pre_ll) {
	bool converged = false;
	double thres_EM = config.getRealPara("thres_EM");
	if(ll-pre_ll < -0.1) {
		//if(config.isVerbose()) {
			cerr << "******likelihood decreased from " << pre_ll << " to " << ll << endl;
		//}
	}
	double ZERO_FINAL = config.getRealPara("ZERO_FINAL");
	double delta_ll = fabs(ll-pre_ll);
	double avg_ll = (fabs(ll)+fabs(pre_ll)+ZERO_FINAL)/2;
	if(delta_ll/avg_ll < thres_EM) {
		converged = true;
	}
	return converged;
}

double Model::estimateParas() {
	double pre_loglik = -1.0e15;
	bool converged = 0;
	int i, num_iter = 1;
	
	int max_iter = config.getIntPara("em_max_iter");
	int verbose = config.getIntPara("verbose");
	while(num_iter <= max_iter && !converged) {
		compute_ess();
		converged = em_converged(loglik, pre_loglik);
		
		if(paras.indicator.get(0,0)) {
			norm_trans(exp_num_visits, 0);
			hmm.setPI(exp_num_visits);
		}
		if(paras.indicator.get(0,1)) {
			double clamp_thres = config.getRealPara("clamp_thres");
			norm_trans(exp_num_trans, clamp_thres);
			hmm.setTransMat(exp_num_trans);
		}
		if(paras.indicator.get(0,2)) {
			paras.lambda = n_paras.lambda;
		}
		if(paras.indicator.get(0,3)) {
			paras.mu = n_paras.mu;
		}
		if(paras.indicator.get(0,4)) {
			paras.eta = n_paras.eta;
		}
		if(paras.indicator.get(0,5)) {
			paras.p = n_paras.p;
		}
		
		/*
		if(verbose) {
			printParas();
			printf("iteration %d, loglik = %.2f\n", num_iter, loglik);
		}
		*/
		
		num_iter++;
		pre_loglik = loglik;
	}
	return loglik;
}

void Model::printParas() {
	int i;
	cerr << "lambda: " << paras.lambda << endl;
	cerr << "mu: " << paras.mu.get(0,0);
	for(i = 1; i < paras.mu.getCOLS(); i++) {
		cerr << ", " << paras.mu.get(0,i);
	}
	cerr << endl;
	cerr << "eta: " << paras.eta.get(0,0);
	for(i = 1; i < paras.eta.getCOLS(); i++) {
		cerr << ", " << paras.eta.get(0,i);
	}
	cerr << endl;
	cerr << "p: " << paras.p.get(0,0);
	for(i = 1; i < paras.p.getCOLS(); i++) {
		cerr << ", " << paras.p.get(0,i);
	}
	cerr << endl;
}

void Model::compute_ess() {
	int i, j, k;
	vector<int>& chromosomes = genomedata.getChroms();
	int num_chr = chromosomes.size();
	
	//initialize
	loglik = 0;
	gamma_sep.clear();
	condi_probs_sep.clear();
	condi_probs_fluct_sep.clear();
	exp_num_visits.set(0);
	exp_num_trans.set(0);
	
	//FwdBack per chr
	vector<int*> threadParas;
	for(i = 0; i < num_chr; i++) {
		int* tparas = new int[2];
		tparas[0] = modelIndx;
		tparas[1] = chromosomes[i];
		threadpool->pool_add_work(&Model::fwdBackPerChr, tparas, i);
		threadParas.push_back(tparas);
	}
	threadpool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
	threadParas.clear();
	
	//update parameters
	n_paras = paras;
	k = paras.mu.getCOLS();
	//update mu
	for(i = 0; i < k; i++) {
		int* tparas = new int[2];
		tparas[0] = modelIndx;
		tparas[1] = i;
		threadpool->pool_add_work(&Model::update_mu, tparas, i);
		threadParas.push_back(tparas);
	}
	//update lambda
	threadpool->pool_add_work(&Model::update_lambda, &modelIndx, 1);
	threadpool->wait();
	
	//update eta
	for(i = 0; i < k; i++) {
		threadpool->pool_add_work(&Model::update_eta, threadParas[i], i);
	}
	
	//update p
	k = paras.p.getCOLS();
	for(i = 0; i < k; i++) {
		int* tparas = new int[2];
		tparas[0] = modelIndx;
		tparas[1] = i;
		threadpool->pool_add_work(&Model::update_p, tparas, i);
		threadParas.push_back(tparas);
	}
	threadpool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
}

void* Model::fwdBackPerChr(const void *arg) {
	int *tparas = (int*) arg;
	int mindx = tparas[0];
	int chr = tparas[1];
	Matrix<double> &chrData = genomedata.getChrData(chr);
	double *p1 = chrData.getEntrance();
	int T = chrData.getCOLS();
	Model &model = cnacaller.getModel(mindx);
	HMM &hmm = model.getHMM();
	
	Matrix<double> tmp;
	vector<Matrix<double> > results(3, tmp);
	vector<Matrix<double> > fbParas(4, tmp);
	model.getObslik(T, p1+T, p1+3*T, results);
	double ll = hmm.FwdBack(results[0], fbParas, 2);
	
	pthread_mutex_t& pm = model.getPM();
	pthread_mutex_lock(&pm);
	model.loglik += ll;
	int num_state = hmm.getStateCount();
	for(int k = 0; k < num_state; k++) {
		model.exp_num_visits.set(0, k, model.exp_num_visits.get(0, k)+fbParas[2].get(k, 0));
	}
	model.exp_num_trans += fbParas[3];
	model.gamma_sep[chr] = fbParas[2];
	model.condi_probs_sep[chr] = results[1];
	model.condi_probs_fluct_sep[chr] = results[2];
	pthread_mutex_unlock(&pm);
}

void* Model::update_mu(const void *arg) {
	int *tparas = (int*) arg;
	int mindx = tparas[0];
	int indx = tparas[1];
	int i, j, k, t;
	Model &model = cnacaller.getModel(mindx);
	vector<int>& chromosomes = genomedata.getChroms();
	int chr, num_chr = chromosomes.size();
	
	double mu = model.n_paras.mu.get(0, indx);
	double eta = model.n_paras.eta.get(0, indx);
	double mu_u = mu;
	double mu_tol = 0.005;
	
	Matrix<double> &mafs = model.hmm.getMAFs();
	double sigma = 0.5;
	
	if(fabs(mafs.get(0,indx)-0.5) > 1.0e-5) {
		int iter = 0;
		double ELL_D_1, ELL_D_2;
		while(1) {
			ELL_D_1 = 0;
			ELL_D_2 = 0;
			for(i = 0; i < num_chr; i++) {
				chr = chromosomes[i];
				if(chr == 23 || chr == 24) {
					continue;
				}
				Matrix<double> &chrData = genomedata.getChrData(chr);
				int T = chrData.getCOLS();
				double *p1 = chrData.getEntrance();
				double *p2 = model.gamma_sep[chr].getEntrance();
				double *p3 = model.condi_probs_sep[chr].getEntrance();
				
				double sum1 = 0, sum2 = 0;
				for(t = 0; t < T; t++) {
					sum1 += p2[indx*T+t]*p3[indx*T+t]*(psi(0,p1[T+t]+mu/eta)-psi(0,p1[2*T+t]-p1[T+t]+(1-mu)/eta)-psi(0,mu/eta)+psi(0,(1-mu)/eta));
					sum2 += p2[indx*T+t]*p3[indx*T+t]*(psi(1,p1[T+t]+mu/eta)+psi(1,p1[2*T+t]-p1[T+t]+(1-mu)/eta)-psi(1,mu/eta)-psi(1,(1-mu)/eta));
				}
				ELL_D_1 += sum1/eta;
				ELL_D_2 += sum2/(eta*eta);
			}
			ELL_D_1 -= (mu-mafs.get(0,indx))/(sigma*sigma);
			ELL_D_2 -= 1/(sigma*sigma);
			
			double mu_adj = -ELL_D_1/ELL_D_2;
			double c_t;
			if(fabs(mafs.get(0,indx)-1.0) < 1.0e-5) {
				c_t = 0.17;
			}
			else {
				c_t = 0.05;
			}
			double c_p = 1;
			if(mu_adj < 0) {
				double temp = (max(0.5,mafs.get(0,indx)-c_t)-mu)/mu_adj;
				if(temp < 1) {
					c_p = temp;
				}
			}
			else if(mu_adj > 0) {
				double temp = (min(0.97,mafs.get(0,indx)+c_t)-mu)/mu_adj;
				if(temp < 1) {
					c_p = temp;
				}
			}
			
			mu_u = mu+c_p*mu_adj;
			if(isnan(mu_u)) {
				mu_u = mu;
			}
			
			iter++;
			
			if(fabs(mu_u-mu) < mu_tol || iter > nr_max_iter) {
				break;
			}
			else {
				mu = mu_u;
			}
		}
	}
	model.n_paras.mu.set(0,indx,mu_u);
}

void* Model::update_lambda(const void *arg) {
	int mindx = *((int*) arg);
	int i, j, k, t;
	Model &model = cnacaller.getModel(mindx);
	vector<int>& chromosomes = genomedata.getChroms();
	int chr, num_chr = chromosomes.size();
	
	double lambda = model.n_paras.lambda;
	double lambda_u = lambda;
	double lambda_tol = 0.5;
	
	int num_state = model.hmm.getStateCount();
	Matrix<double> &cns = model.hmm.getCNs();
	Matrix<double> p_temp(1, num_state);
	Matrix<double> tmp1(1, num_state);
	for(i = 0; i < num_state; i++) {
		k = round(cns.get(0, i));
		double v = model.n_paras.p.get(0,k);
		p_temp.set(0, i, v);
		tmp1.set(0, i, (1-v)/v);
	}
	
	int iter = 0;
	double ELL_D_1, ELL_D_2;
	while(1) {
		ELL_D_1 = 0;
		ELL_D_2 = 0;
		Matrix<double> lambda_c = cns*(lambda/2);
		for(i = 0; i < num_chr; i++) {
			chr = chromosomes[i];
			if(chr == 23 || chr == 24) {
				continue;
			}
			Matrix<double> &chrData = genomedata.getChrData(chr);
			int T = chrData.getCOLS();
			double *p1 = chrData.getEntrance();
			double *p2 = model.gamma_sep[chr].getEntrance();
			double *p3 = model.condi_probs_fluct_sep[chr].getEntrance();
			
			double sum1 = 0, sum2 = 0;
			for(j = 0; j < num_state; j++) {
				for(t = 0; t < T; t++) {
					double post_probs = p2[j*T+t]*(1-p3[j*T+t]);
					double v1 = cns.get(0,j)*tmp1.get(0,j)/2;
					double v2 = psi(0,p1[3*T+t]+lambda_c.get(0,j)*tmp1.get(0,j))+log(1-p_temp.get(0,j))-psi(0,lambda_c.get(0,j)*tmp1.get(0,j));
					double v3 = psi(1,p1[3*T+t]+lambda_c.get(0,j)*tmp1.get(0,j))-psi(1,lambda_c.get(0,j)*tmp1.get(0,j));
					sum1 += post_probs*v1*v2;
					sum2 += post_probs*v1*v1*v3;
				}
			}
			ELL_D_1 += sum1;
			ELL_D_2 += sum2;
		}
		
		double lambda_adj = -ELL_D_1/ELL_D_2;
		double c_p = 1;
		if(lambda_adj < 0) {
			double temp = (model.lambda_limit.get(0,0)-lambda)/lambda_adj;
			if(temp < 1) {
				c_p = temp;
			}
		}
		else if(lambda_adj > 0) {
			double temp = (model.lambda_limit.get(0,1)-lambda)/lambda_adj;
			if(temp < 1) {
				c_p = temp;
			}
		}
		
		lambda_u = lambda+c_p*lambda_adj;
		if(isnan(lambda_u)) {
			lambda_u = lambda;
		}
		//cerr << "lambda: " << lambda_u << endl;
		
		iter++;
		
		if(fabs(lambda_u-lambda) < lambda_tol || iter > nr_max_iter) {
			break;
		}
		else {
			lambda = lambda_u;
		}
	}
	model.n_paras.lambda = lambda_u;
}

void* Model::update_eta(const void *arg) {
	int *tparas = (int*) arg;
	int mindx = tparas[0];
	int indx = tparas[1];
	int i, j, k, t;
	Model &model = cnacaller.getModel(mindx);
	vector<int>& chromosomes = genomedata.getChroms();
	int chr, num_chr = chromosomes.size();
	
	double mu = model.n_paras.mu.get(0, indx);
	double eta = model.n_paras.eta.get(0, indx);
	double eta_u = eta;
	double eta_tol = 0.01;
	
	int iter = 0;
	double ELL_D_1, ELL_D_2;
	while(1) {
		ELL_D_1 = 0;
		ELL_D_2 = 0;
		for(i = 0; i < num_chr; i++) {
			chr = chromosomes[i];
			if(chr == 23 || chr == 24) {
				continue;
			}
			Matrix<double> &chrData = genomedata.getChrData(chr);
			int T = chrData.getCOLS();
			double *p1 = chrData.getEntrance();
			double *p2 = model.gamma_sep[chr].getEntrance();
			double *p3 = model.condi_probs_sep[chr].getEntrance();
			
			double sum1 = 0, sum2 = 0;
			for(t = 0; t < T; t++) {
				double post_probs = p2[indx*T+t]*p3[indx*T+t];
				double v1 = -mu*(psi(0,p1[T+t]+mu/eta)-psi(0,p1[2*T+t]+1/eta));
				double v2 = -(1-mu)*(psi(0,p1[2*T+t]-p1[T+t]+(1-mu)/eta)-psi(0,p1[2*T+t]+1/eta));
				double v3 = mu*(psi(0,mu/eta)-psi(0,1/eta));
				double v4 = (1-mu)*(psi(0,(1-mu)/eta)-psi(0,1/eta));
				double v1_1 = mu*(mu*psi(1,p1[T+t]+mu/eta)-psi(1,p1[2*T+t]+1/eta));
				double v2_1 = (1-mu)*((1-mu)*psi(1,p1[2*T+t]-p1[T+t]+(1-mu)/eta)-psi(1,p1[2*T+t]+1/eta));
				double v3_1 = -mu*(mu*psi(1,mu/eta)-psi(1,1/eta));
				double v4_1 = -(1-mu)*((1-mu)*psi(1,(1-mu)/eta)-psi(1,1/eta));
				sum1 += post_probs*(v1+v2+v3+v4)/pow(eta,2);
				sum2 += post_probs*(-(v1+v2+v3+v4)*2/pow(eta,3)+(v1_1+v2_1+v3_1+v4_1)/pow(eta,4));
			}
			ELL_D_1 += sum1;
			ELL_D_2 += sum2;
		}
		
		double eta_adj = -ELL_D_1/ELL_D_2;
		double c_p = 1;
		if(eta_adj < 0) {
			double temp = (0.001-eta)/eta_adj;
			if(temp < 1) {
				c_p = temp;
			}
		}
		else if(eta_adj > 0) {
			double temp = (1000-eta)/eta_adj;
			if(temp < 1) {
				c_p = temp;
			}
		}
		
		eta_u = eta+c_p*eta_adj;
		if(isnan(eta_u)) {
			eta_u = eta;
		}
		
		//cerr << "eta: " << eta << endl;
		
		iter++;
		
		if(fabs(eta_u-eta) < eta_tol || iter > nr_max_iter) {
			break;
		}
		else {
			eta = eta_u;
		}
	}
	model.n_paras.eta.set(0,indx,eta_u);
}

void* Model::update_p(const void *arg) {
	int *tparas = (int*) arg;
	int mindx = tparas[0];
	int indx = tparas[1];
	int i, j, k, t;
	Model &model = cnacaller.getModel(mindx);
	vector<int>& chromosomes = genomedata.getChroms();
	int chr, num_chr = chromosomes.size();
	
	double lambda = model.n_paras.lambda;
	double p = model.n_paras.p.get(0, indx);
	double p_u = p;
	double p_tol = 1.0e-5;
	
	int num_state = model.hmm.getStateCount();
	Matrix<double> &cns = model.hmm.getCNs();
	Matrix<double> lambda_c = cns*(lambda/2);
	
	int iter = 0;
	double ELL_D_1, ELL_D_2;
	while(1) {
		ELL_D_1 = 0;
		ELL_D_2 = 0;
		Matrix<double> p_temp(1, num_state);
		Matrix<double> tmp1(1, num_state);
		for(i = 0; i < num_state; i++) {
			k = round(cns.get(0, i));
			if(k == indx) {
				p_temp.set(0, i, p);
				tmp1.set(0, i, lambda_c.get(0,i)*(1-p)/p);
			}
		}
		for(i = 0; i < num_chr; i++) {
			chr = chromosomes[i];
			if(chr == 23 || chr == 24) {
				continue;
			}
			Matrix<double> &chrData = genomedata.getChrData(chr);
			int T = chrData.getCOLS();
			double *p1 = chrData.getEntrance();
			double *p2 = model.gamma_sep[chr].getEntrance();
			double *p3 = model.condi_probs_fluct_sep[chr].getEntrance();
			
			double sum1 = 0, sum2 = 0;
			for(j = 0; j < num_state; j++) {
				k = round(cns.get(0, j));
				if(k != indx) {
					continue;
				}
				for(t = 0; t < T; t++) {
					double post_probs = p2[j*T+t]*(1-p3[j*T+t]);
					double v1 = lambda_c.get(0,j)/pow(p_temp.get(0,j),2)*(psi(0,tmp1.get(0,j))-psi(0,p1[3*T+t]+tmp1.get(0,j))-log(1-p_temp.get(0,j)));
					double v2 = (p1[3*T+t]-lambda_c.get(0,j))/p_temp.get(0,j);
					double v1_1 = 2*lambda_c.get(0,j)/p_temp.get(0,j)*(psi(0,p1[3*T+t]+tmp1.get(0,j))+log(1-p_temp.get(0,j))-psi(0,tmp1.get(0,j)));
					double v2_1 = lambda_c.get(0,j)*(lambda_c.get(0,j)/pow(p_temp.get(0,j),2)*psi(1,p1[3*T+t]+tmp1.get(0,j))+1/(1-p_temp.get(0,j)));
					double v3_1 = lambda_c.get(0,j)*(-lambda_c.get(0,j)/pow(p_temp.get(0,j),2)*psi(1,tmp1.get(0,j)))-p1[3*T+t]+lambda_c.get(0,j);
					sum1 += post_probs*(v1+v2);
					sum2 += post_probs/pow(p_temp.get(0,j),2)*(v1_1+v2_1+v3_1);
				}
			}
			ELL_D_1 += sum1;
			ELL_D_2 += sum2;
		}
		
		double p_adj = -ELL_D_1/ELL_D_2;
		double c_p = 1;
		if(p_adj < 0) {
			double temp = (1.0e-5-p)/p_adj;
			if(temp < 1) {
				c_p = temp;
			}
		}
		else if(p_adj > 0) {
			double temp = (1-1.0e-5-p)/p_adj;
			if(temp < 1) {
				c_p = temp;
			}
		}
		
		p_u = p+c_p*p_adj;
		if(isnan(p_u)) {
			p_u = p;
		}
		
		iter++;
		
		if(fabs(p_u-p) < p_tol || iter > nr_max_iter) {
			break;
		}
		else {
			p = p_u;
		}
	}
	model.n_paras.p.set(0,indx,p_u);
}

void Model::getObslik(Matrix<double> &rd, Matrix<double> &rc, vector<Matrix<double> > &results) {
	int i, j, t, indx;
	double cn, maf;
	int N = hmm.getPI().getCOLS(); //number of states
	int T = rd.getCOLS(); //number of observations
	
	Matrix<double> &cns = hmm.getCNs();
	Matrix<double> &mafs = hmm.getMAFs();
	
	double fluct_prob = 0.005;
	double het_prior = 2.0/3;
	double *p1 = rd.getEntrance();
	double *p2 = rc.getEntrance();
	
	double max_rc = genomedata.getMaxRC();
	Matrix<double> alphas = paras.mu.dotDivide(paras.eta);
	Matrix<double> betas(1, N, false); 
	for(i = 0; i < N; i++) {
		betas.set(0, i, (1-paras.mu.get(0,i))/paras.eta.get(0,i));
	}
	
	results[0].resize(N, T, false); // obslik
	results[1].resize(N, T, false); // condi_probs
	results[2].resize(N, T, false); // condi_probs_fluct
	
	double *p3 = results[0].getEntrance();
	double *p4 = results[1].getEntrance();
	double *p5 = results[2].getEntrance();
	
	Matrix<double> obslik_md_Homo(1, T, false);
	double *p6 = obslik_md_Homo.getEntrance();
	for(t = 0; t < T; t++) {
		p6[t] = binopdf(p1[t], p1[T+t], 0.9999);
	}
	
	double prior_p;
	for(i = 0; i < N; i++) {
		cn = cns.get(0,i);
		maf = mafs.get(0,i);
		if(fabs(maf-0.5) < 1.0e-5) {
			prior_p = paras.normal_prior.get(0, cn/2);
		}
		else {
			prior_p = 1;
		}
		for(t = 0; t < T; t++) {
			double obslik_rc = hmm.getObslik_rc(p2[t], paras.lambda*cn/2, paras.p.get(0,round(cn)));
			double obslik_md_Het = prior_p*hmm.getObslik_md(p1[t], p1[T+t], alphas.get(0,i), betas.get(0,i));
			double obslik_md = (1-het_prior)*p6[t]+het_prior*obslik_md_Het;
			
			p3[i*T+t] = (1-fluct_prob)*obslik_rc*obslik_md+fluct_prob/((p1[T+t]+1)*max_rc);
			p4[i*T+t] = (1-fluct_prob)*obslik_rc*het_prior*obslik_md_Het/p3[i*T+t];
			p5[i*T+t] = (fluct_prob/((p1[T+t]+1)*max_rc))/p3[i*T+t];
		}
	}
}

void Model::getObslik(int T, double *rd, double *rc, vector<Matrix<double> > &results) {
	int i, j, t, indx;
	double cn, maf;
	int N = hmm.getPI().getCOLS(); //number of states
	
	Matrix<double> &cns = hmm.getCNs();
	Matrix<double> &mafs = hmm.getMAFs();
	
	double fluct_prob = 0.005;
	double het_prior = 2.0/3;
	
	double max_rc = genomedata.getMaxRC();
	Matrix<double> alphas = paras.mu.dotDivide(paras.eta);
	Matrix<double> betas(1, N, false); 
	for(i = 0; i < N; i++) {
		betas.set(0, i, (1-paras.mu.get(0,i))/paras.eta.get(0,i));
	}
	
	results[0].resize(N, T, false); // obslik
	results[1].resize(N, T, false); // condi_probs
	results[2].resize(N, T, false); // condi_probs_fluct
	double *p1 = results[0].getEntrance();
	double *p2 = results[1].getEntrance();
	double *p3 = results[2].getEntrance();
	
	Matrix<double> obslik_md_Homo(1, T, false);
	double *p4 = obslik_md_Homo.getEntrance();
	for(t = 0; t < T; t++) {
		p4[t] = binopdf(rd[t], rd[T+t], 0.9999);
	}
	
	double prior_p;
	for(i = 0; i < N; i++) {
		cn = cns.get(0,i);
		maf = mafs.get(0,i);
		if(fabs(maf-0.5) < 1.0e-5) {
			prior_p = paras.normal_prior.get(0, cn/2);
		}
		else {
			prior_p = 1;
		}
		for(t = 0; t < T; t++) {
			double obslik_rc = hmm.getObslik_rc(rc[t], paras.lambda*cn/2, paras.p[round(cn)]);
			double obslik_md_Het = prior_p*hmm.getObslik_md(rd[t], rd[T+t], alphas[i], betas[i]);
			double obslik_md = (1-het_prior)*p4[t]+het_prior*obslik_md_Het;
			
			p1[i*T+t] = (1-fluct_prob)*obslik_rc*obslik_md+fluct_prob/((rd[T+t]+1)*max_rc);
			p2[i*T+t] = (1-fluct_prob)*obslik_rc*het_prior*obslik_md_Het/p1[i*T+t];
			p3[i*T+t] = (fluct_prob/((rd[T+t]+1)*max_rc))/p1[i*T+t];
		}
	}
}

void Model::process_results(Matrix<double> &p_states, double &acn, map<int, Matrix<double> > &segments_sep) {
	int i, j, k;
	int chr;
	vector<int>& chroms = genomedata.getChroms();
	
	int sindx, eindx;
	int num_state = hmm.getStateCount();
	Matrix<double> exp_num_states(num_state, 1, 0);
	Matrix<int> stateseq;
	double cn, weight_len = 0, total_len = 0;
	segments_sep.clear();
	int num_snp = 0;
	for(i = 0; i < chroms.size(); i++) {
		chr = chroms[i];
		Matrix<double> &chrData = genomedata.getChrData(chr);
		double *p1 = chrData.getEntrance();
		Matrix<double> &post_probs = gamma_sep[chr];
		Matrix<double> tmp = post_probs.sumCols();
		exp_num_states += tmp;
		num_snp += chrData.getCOLS();
		
		post_probs.max(1, stateseq);
		getSegments(stateseq, segments_sep[chr]);
		int num_seg = segments_sep[chr].getROWS();
		for(j = 0; j < num_seg; j++) {
			sindx = segments_sep[chr].get(j,1);
			eindx = segments_sep[chr].get(j,2);
			cn = segments_sep[chr].get(j,3);
			weight_len += cn*(p1[eindx]-p1[sindx]+1);
			total_len += p1[eindx]-p1[sindx]+1;
		}
	}
	
	p_states = exp_num_states/num_snp;
	acn = weight_len/total_len;
}

void Model::predictStates() {
	int i;
	vector<int>& chroms = genomedata.getChroms();
	vector<int*> threadParas;
	for(i = 0; i < chroms.size(); i++) {
		int* tparas = new int[2];
		tparas[0] = modelIndx;
		tparas[1] = chroms[i];
		threadpool->pool_add_work(&Model::predictChrStates, tparas, i);
		threadParas.push_back(tparas);
	}
	threadpool->wait();
	for(i = 0; i < threadParas.size(); i++) {
		delete[] threadParas[i];
	}
}

void* Model::predictChrStates(const void *arg) {
	int *tparas = (int*) arg;
	int mindx = tparas[0];
	int chr = tparas[1];

	Matrix<double> &chrData = genomedata.getChrData(chr);
	double *p1 = chrData.getEntrance();
	int T = chrData.getCOLS();
	Model &model = cnacaller.getModel(mindx);
	HMM &hmm = model.getHMM();
	
	Matrix<double> tmp;
	vector<Matrix<double> > results(3, tmp);
	model.getObslik(T, p1+T, p1+3*T, results);
	
	Matrix<int> stateseq;
	Matrix<double> segments;
	hmm.Viterbi(results[0], stateseq); // viterbi algorithm
	model.getSegments(stateseq, segments);
	model.correctStates(chr, segments);
	model.calculateScores(chr, segments);
	genomedata.setSegments(chr, segments);
}

void Model::getSegments(Matrix<int> &stateseq, Matrix<double> &segments) {
	int i, k = 0;
	int pre_state = -1;
	int s_indx = -1;
	int T = stateseq.getCOLS();
	int *p1 = stateseq.getEntrance();
	segments.resize(T, 6, false);
	double *p2 = segments.getEntrance();
	
	double *p3 = hmm.getCNs().getEntrance();
	double *p4 = hmm.getMAFs().getEntrance();
	
	for(i = 0; i < T; i++) {
		if(pre_state == -1) {
			pre_state = p1[i];
			s_indx = i;
		}
		else if(p1[i] != pre_state) {
			p2[k*6] = pre_state;
			p2[k*6+1] = s_indx;
			p2[k*6+2] = i-1;
			p2[k*6+3] = p3[pre_state];
			p2[k*6+4] = round(p3[pre_state]*p4[pre_state]);
			k++;
			pre_state = p1[i];
			s_indx = i;
		}
	}
	
	p2[k*6] = pre_state;
	p2[k*6+1] = s_indx;
	p2[k*6+2] = T-1;
	p2[k*6+3] = p3[pre_state];
	p2[k*6+4] = round(p3[pre_state]*p4[pre_state]);
	
	segments.resize(k+1, 6, true);
}

void Model::correctStates(int chr, Matrix<double> &segments) {
	int i, j, k, t;
	int num_segs = segments.getROWS();
	
	Matrix<double> chrData = genomedata.getChrData(chr);
	double *p1 = chrData.getEntrance();
	int T = chrData.getCOLS();
	
	double *p2 = segments.getEntrance();
	
	int max_cn = config.getIntPara("maxCN");
	Matrix<double> &cns = hmm.getCNs();
	Matrix<double> &mafs = hmm.getMAFs();
	
	int state, sindx, eindx;
	double cn, cn_p, mcn;
	double prob;
	double homo_th = config.getRealPara("homo_th");
	double maf, median_rc, median_maf;
	for(i = 0; i < num_segs; i++) {
		state = p2[i*6];
		sindx = p2[i*6+1];
		eindx = p2[i*6+2];
		cn = p2[i*6+3];
		
		median_rc = median(&p1[3*T+sindx], eindx-sindx+1);
		cn_p = round(2*median_rc/paras.lambda);
		if(cn_p <= max_cn+0.01 && fabs(cn-cn_p) < 0.001) { // state is included in the HMM and accurately detected
			continue;
		}
		p2[i*6] = -1;
		p2[i*6+3] = cn_p;
		
		vector<double> mafs_het;
		for(t = sindx; t <= eindx; t++) {
			prob = binopdf(p1[T+t], p1[2*T+t], 0.01);
			if(prob >= homo_th) {
				continue;
			}
			maf = p1[T+t]/p1[2*T+t];
			mafs_het.push_back(maf);
		}
		if(mafs_het.empty()) {
			p2[i*6+4] = cn_p;
		}
		else {
			median_maf = median(mafs_het);
			double min_dist = 1;
			int minIndx;
			for(j = ceil((cn_p-0.5)/2); j <= cn_p; j++) {
				double dist = fabs(j/cn_p-median_maf);
				if(dist < min_dist) {
					min_dist = dist;
					minIndx = j;
				}
			}
			p2[i*6+4] = minIndx;
		}
		if(cn_p < max_cn+0.01) {
			maf = p2[i*6+4]/cn_p;
			for(j = 0; j < cns.getCOLS(); j++) {
				if(fabs(cn_p-cns[j])<0.02 && fabs(maf-mafs[j])<0.001) {
					p2[i*6] = j;
					break;
				}
			}
		}
	}
}

void Model::calculateScores(int chr, Matrix<double> &segments) {
	int i, j, t;
	int num_segs = segments.getROWS();
	
	Matrix<double> &chrData = genomedata.getChrData(chr);
	double *p1 = chrData.getEntrance();
	int T = chrData.getCOLS();
	
	double *p2 = segments.getEntrance();
	
	Matrix<double> &cns = hmm.getCNs();
	Matrix<double> &mafs = hmm.getMAFs();
	
	int state, sindx, eindx;
	double cn, mcn;
	double prob;
	double homo_th = config.getRealPara("homo_th");
	double maf, maf_e, score;
	double lambda_e;
	for(i = 0; i < num_segs; i++) {
		state = p2[i*6];
		sindx = p2[i*6+1];
		eindx = p2[i*6+2];
		cn = p2[i*6+3];
		mcn = p2[i*6+4];
		
		if(state < 0) {
			maf_e = mcn/cn;
		}
		else {
			maf_e = paras.mu.get(0, state);
		}
		lambda_e = cn*paras.lambda/2;
		vector<double> scores;
		for(t = sindx; t <= eindx; t++) {
			prob = binopdf(p1[T+t], p1[2*T+t], 0.01);
			if(prob >= homo_th) {
				continue;
			}
			maf = p1[T+t]/p1[2*T+t];
			score = exp(-fabs(maf-maf_e)/maf_e);
			score = 100*(score+exp(-fabs(p1[3*T+t]-lambda_e)/lambda_e))/2;
			scores.push_back(score);
		}
		if(scores.empty()) {
			p2[i*6+5] = -1;
		}
		else {
			p2[i*6+5] = median(scores);
		}
	}
}

