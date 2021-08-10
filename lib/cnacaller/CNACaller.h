// ***************************************************************************
// CNACaller.h (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _CNACALLER_H
#define _CNACALLER_H

#include <pthread.h>

#include "Model.h"
#include "Matrix.h"

class CNACaller {
	private:
		Model *models;
		int num_model;
		int best_indx;
		bool noSolutionFlag;
		
		Matrix<double> LL_all;
		Matrix<double> acn_all;
		Matrix<double> p_states_all;
		Matrix<double> p_total_del;
		Matrix<int> valid_indic;
		
		pthread_mutex_t pm;
		
		void init();
		void getBestModel();
		void predict();
		void plotResults();
		void plotResults_new();
		
		void inferModelParas(int mindx);
		
		void plotChrResults(int chr);
		static void* plotChrResults(const void *arg);
		
	public:
		Model& getModel(int indx) {return models[indx];}
		void callCNAs();
		
};


#endif

