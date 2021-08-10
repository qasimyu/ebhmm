// ***************************************************************************
// GenomeData.cpp (c) 2019 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <cassert>

#include "GenomeData.h"
#include "Config.h"
#include "MyDefine.h"
#include "Fasta.h"
#include "split.h"

#include "api/BamReader.h"
#include "api/BamWriter.h"

extern "C" {
	#include "common.h"
	#include "bbiFile.h"
	#include "bigWig.h"
}

using namespace std;
using namespace BamTools;

void GenomeData::LoadAndProcess() {
	/*** parameters for filtering SNPs ***/
	double min_td = config.getRealPara("min_td");
	double max_td = config.getRealPara("max_td");
	double min_gc = config.getRealPara("min_gc");
	double min_map = config.getRealPara("min_map");
	double max_map = config.getRealPara("max_map");
	double rc_th[2] = {0.1, 99.5};
	
	/*** fetch sample name ***/
	string inputFile = config.getStringPara("input");
	int indx1 = inputFile.find_last_of('/');
	if(indx1 == string::npos) {
		indx1 = -1;
	}
	int indx2 = inputFile.find_last_of('.');
	if(indx2 == string::npos) {
		indx2 = inputFile.length();
	}
	sampleName = inputFile.substr(indx1+1, indx2-indx1-1);
	
	/*** load data ***/
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open input file " << inputFile << endl;
		exit(-1);
	}
	fgets(buf, 500, fp);
	fclose(fp);
	int num_snp = atoi(buf)-2;
	
	Matrix<double> data_all(7, num_snp, false); // chr, pos, bd, td, rc, gc ,map
	if(!(fp = fopen(inputFile.c_str(),"r"))){
		cerr << "Error: cannot open input file " << inputFile << endl;
		exit(-1);
	}
	
	char* elems[10];
	int elemnum;
	long line_num = 2;
	int i = 0, j, k;
	fgets(buf, 500, fp);
	fgets(buf, 500, fp); //skip first two lines
    while(fgets(buf, 500, fp)) {
		line_num++;
		elemnum = split(buf, '\t', elems);
		if(elemnum < 7){
			cerr << "Warning: malformed input file " << inputFile <<
                    ", there should be 7 fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(j = 0; j < 7; j++) {
			data_all.set(j, i, atof(elems[j]));
		}
		i++;
	}
	fclose(fp);
	assert(i == num_snp);
	
	cerr << "Total " << num_snp << " SNPs are loaded from file " << inputFile << endl;
	if(num_snp < 100000) {
		cerr << "Warning: the number of SNPs loaded is too small, check the format of data file and whether data are completely loaded!" << endl;
	}
	
	/*** initial filtering ***/
	double *p1 = data_all.getEntrance();
	vector<int> indxs;
	double rc_low_p = prctile(p1+4*num_snp, num_snp, rc_th[0]);
	double rc_high_p = prctile(p1+4*num_snp, num_snp, rc_th[1]);
	for(i = 0; i < num_snp; i++) {
		bool flag = p1[5*num_snp+i] > min_gc; // gc filtering
		flag = flag && p1[6*num_snp+i] > min_map && p1[6*num_snp+i] < max_map; // map filtering
		//flag = flag && p1[4*num_snp+i] > rc_low_p && p1[4*num_snp+i] < rc_high_p; // rc filtering
		flag = flag && p1[4*num_snp+i] < rc_high_p; // rc filtering
		flag = flag && p1[3*num_snp+i] >= min_td && p1[3*num_snp+i] <= max_td; // td filtering
		if(flag) {
			indxs.push_back(i);
		}
	}
	data_all = data_all.Cols(indxs);
	indxs.clear();
	num_snp = data_all.getCOLS();
	p1 = data_all.getEntrance();
	
	/*** correct GC-content and mappability bias ***/
	corrRC(p1+4*num_snp, p1+6*num_snp, num_snp); // Map-Correction
	corrRC(p1+4*num_snp, p1+5*num_snp, num_snp);// GC-Correction
	
	/*** filter SNPs with extremely high rc after correction ***/
	rc_high_p = prctile(p1+4*num_snp, num_snp, 99.99);
	for(i = 0; i < num_snp; i++) {
		if(p1[4*num_snp+i] < rc_high_p) {
			indxs.push_back(i);
		}
	}
	data_all = data_all.Cols(indxs);
	indxs.clear();
	num_snp = data_all.getCOLS();
	p1 = data_all.getEntrance();
	
	/*** remove bias in read depth data ***/
	// Quantile normalization of BAF
	//preprocessRD(p1+2*num_snp, p1+3*num_snp, num_snp);
	
	/*** filter homozygous positions, #het/#homo=2/1 ***/
	filterSNPs(data_all);
	num_snp = data_all.getCOLS();
	p1 = data_all.getEntrance();
	
	double *tmp = new double[num_snp];
	memcpy(tmp, p1+4*num_snp, sizeof(double)*num_snp);
	rc_median = median(tmp, num_snp);
	max_rc = tmp[num_snp-1];
	delete[] tmp;
	double clamp_thres = 1-1.0e-4/num_snp;
	config.setRealPara("clamp_thres", clamp_thres);
	
	/*** use at least 100000 data points for screening ***/
	stepsize_ds = max(1, num_snp/100000);
	
	/*** group by chromosomes ***/
	splitDataByChr(data_all);
	
	/*** save processed data ***/
	saveData(data_all);
}

void* GenomeData::smoothRCByChr(const void *arg) {
	int chr = *((int*) arg);
	
	int T = genomedata.data_sep[chr].getCOLS();
	double *rc_p = genomedata.data_sep[chr].getEntrance()+3*T;
	int half_w_size = 100;
	int i, j, k;
	int indx;
	double *tmp = new double[2*half_w_size];
	double *rc_n_p = new double[T];
	for(i = 0; i < T; i++) {
		indx = max(0, i-half_w_size);
		j = 0;
		while(indx < i) {
			tmp[j++] = rc_p[indx++];
		}
		indx++;
		while(indx < T) {
			tmp[j++] = rc_p[indx++];
			if(j == half_w_size*2) {
				break;
			}
		}
		double rc = rc_p[i];
		double m_value = median(tmp, j);
		if(fabs(rc-m_value)/m_value > 1) {
			rc_n_p[i] = m_value;
		}
		else {
			rc_n_p[i] = rc;
		}
	}
	
	for(i = 0; i < T; i++) {
		rc_p[i] = rc_n_p[i];
	}
	delete[] tmp;
	delete[] rc_n_p;
}

/*** remove distribution bias of allele read depth using quantile normalization ***/
void GenomeData::preprocessRD(double *bd_p, double *td_p, int T) {
	int i, j, k;
	double *a_fre = new double[T];
	double *b_fre = new double[T];
	double *tmp1 = new double[T];
	double *tmp2 = new double[T];
	for(i = 0, k = 0; i < T; i++) {
		b_fre[i] = bd_p[i]/td_p[i];
		a_fre[i] = 1-b_fre[i];
		if(b_fre[i] > 0.1 && b_fre[i] < 0.9) {
			tmp1[k] = a_fre[i];
			tmp2[k] = b_fre[i];
			k++;
		}
	}
	QuickSort(tmp1, k, false);
	QuickSort(tmp2, k, false);
	double sum = 0;
	for(i = 0; i < k; i++) {
		sum += fabs(tmp1[i]-tmp2[i]);
	}
	/*** average distance before tQN ***/
	double dist1 = sum/k;
	
	int *indx1 = QuickSort(a_fre, T, true);
	int *indx2 = QuickSort(b_fre, T, true);
	double m;
	for(i = 0; i < T; i++) {
		m = (a_fre[i]+b_fre[i])/2;
		tmp1[indx1[i]] = m;
		tmp2[indx2[i]] = m;
		b_fre[i] = bd_p[i]/td_p[i];
		a_fre[i] = 1-b_fre[i];
	}
	delete[] indx1;
	delete[] indx2;
	
	/*** find best threshold value ***/
	double th, eps = config.getRealPara("ZERO_FINAL");
	double *tmp3 = new double[T];
	double *tmp4 = new double[T];
	double a, b, baf;
	double best_th, min_dist = dist1;
	for(th = 0.5; th <= 1.5; th += 0.05) {
		for(i = 0, k = 0; i < T; i++) {
			if(tmp1[i]/(a_fre[i]+eps) > th) {
				a = a_fre[i]*th;
			}
			else {
				a = tmp1[i];
			}
			if(tmp2[i]/(b_fre[i]+eps) > th) {
				b = b_fre[i]*th;
			}
			else {
				b = tmp2[i];
			}
			baf = b/(a+b);
			if(baf > 0.1 && baf < 0.9) {
				tmp3[k] = 1-baf;
				tmp4[k] = baf;
				k++;
			}
		}
		QuickSort(tmp3, k, false);
		QuickSort(tmp4, k, false);
		sum = 0;
		for(i = 0; i < k; i++) {
			sum += fabs(tmp3[i]-tmp4[i]);
		}
		sum /= k;
		if(sum < min_dist) {
			min_dist = sum;
			best_th = th;
		}
	}
	delete[] tmp3;
	delete[] tmp4;
	
	if(dist1 > min_dist) {
		for(i = 0; i < T; i++) {
			if(tmp1[i]/(a_fre[i]+eps) > best_th) {
				tmp1[i] = a_fre[i]*best_th;
			}
			if(tmp2[i]/(b_fre[i]+eps) > best_th) {
				tmp2[i] = b_fre[i]*best_th;
			}
			baf = tmp2[i]/(tmp1[i]+tmp2[i]);
			if(baf > 1e-5 && 1-baf > 1e-5) {
				int d = round((td_p[i]*baf-bd_p[i])/(1-baf));
				bd_p[i] += d;
				td_p[i] += d;
			}
		}
	}
	delete[] a_fre;
	delete[] b_fre;
	delete[] tmp1;
	delete[] tmp2;
}

void GenomeData::filterSNPs(Matrix<double> &data_all) {
	int i, j, k;
	double *p1 = data_all.getEntrance();
	int num_snp = data_all.getCOLS();
	int *status = new int[num_snp];
	double homo_th = config.getRealPara("homo_th");
	double a, b;
	int hetCount = 0;
	for(i = 0; i < num_snp; i++) {
		a = binopdf(p1[2*num_snp+i], p1[3*num_snp+i], 0.01);
		b = binopdf(p1[3*num_snp+i]-p1[2*num_snp+i], p1[3*num_snp+i], 0.01);
		if(a < homo_th && b < homo_th) {
			status[i] = 1;
			hetCount++;
		}
		else {
			status[i] = 0;
		}
	}
	int step = max(1, (int) floor((num_snp-hetCount)/(0.5*hetCount))); // #het/#homo=2/1
	
	vector<int> indxs;
	k = -1;
	for(i = 0; i < num_snp; i++) {
		if(status[i] == 1) {
			indxs.push_back(i);
		}
		else {
			k = (k+1)%step;
			if(k == 0) {
				indxs.push_back(i);
			}
		}
	}
	delete[] status;
	
	data_all = data_all.Cols(indxs);
	
}

void GenomeData::splitDataByChr(Matrix<double> &data_all) {
	int i, j, k;
	int T;
	double *p1 = data_all.getEntrance();
	int num_snp = data_all.getCOLS();
	
	int chr, curChr = -1;
	vector<int> sindxs, eindxs;
	int indx;
	chromosomes.clear();
	for(i = 0; i < num_snp; i++) {
		chr = p1[i];
		if(curChr == -1) {
			chromosomes.push_back(chr);
			sindxs.push_back(i);
			curChr = chr;
		}
		else if(curChr != chr) {
			eindxs.push_back(i-1);
			chromosomes.push_back(chr);
			sindxs.push_back(i);
			curChr = chr;
		}
	}
	eindxs.push_back(i-1);
	
	int tmp[] = {1, 2, 3, 4};
	vector<int> rowIndxs(tmp, tmp+sizeof(tmp)/sizeof(int));
	for(i = 0; i < chromosomes.size(); i++) {
		chr = chromosomes[i];
		data_all.subMat(rowIndxs, sindxs[i], eindxs[i], data_sep[chr]);
		threadpool->pool_add_work(&GenomeData::smoothRCByChr, &chromosomes[i], i);
		T = data_sep[chr].getCOLS();
		double *p2 = data_sep[chr].getEntrance();
		for(j = 0; j < T; j++) {
			double baf = p2[T+j]/p2[2*T+j];
			if(0.5-baf > 1e-5) {
				p2[T+j] = p2[2*T+j]-p2[T+j];
			}
		}	
	}
	threadpool->wait();
	
	for(i = 0; i < chromosomes.size(); i++) {
		chr = chromosomes[i];
		T = data_sep[chr].getCOLS();
		T = (T-1)/stepsize_ds+1;
		data_sep_ds[chr].resize(data_sep[chr].getROWS(), T, false);
		for(j = 0; j < data_sep[chr].getROWS(); j++) {
			for(k = 0; k < T; k++) {
				data_sep_ds[chr].set(j, k, data_sep[chr].get(j, k*stepsize_ds));
			}	
		}
		
	}
}

void GenomeData::corrRC(double *rcs_p, double *factors, int T) {
	double *p1 = rcs_p;
	double *p2 = factors;
	int i, j, t;
	
	char buf[100];
	map<int, vector<int> > pos_sep;
	map<int, vector<int> >::iterator it;
	for(t = 0; t < T; t++) {
		int key = p2[t]*100;
		pos_sep[key].push_back(t);	
	}
	
	Matrix<double> tmp(1, T, false);
	for(i = 0; i < T; i++) {
		tmp.set(0, i, p1[i]);
	}
	double m_all = median(tmp.getEntrance(), T);
	tmp.clear();
	
	for(it = pos_sep.begin(); it != pos_sep.end(); it++) {
		vector<int> pos = (*it).second;
		Matrix<double> f_ratio(1, pos.size() ,false);
		double *p3 = f_ratio.getEntrance();
		for(j = 0; j < pos.size(); j++) {
			t = pos[j];
			p3[j] = p1[t];
		}
		double m_gc = median(p3, pos.size());
		if(m_gc > 0) {
			for(j = 0; j < pos.size(); j++) {
				t = pos[j];
				p1[t] = round(p1[t]*m_all/m_gc);
			}
		}
	}
}

void GenomeData::saveData(Matrix<double> &data_all) {
	int i, j;
	int N = data_all.getROWS();
	int T = data_all.getCOLS();
	
	string outputDir = config.getStringPara("output");
	string outputFile = outputDir+"/"+sampleName+".processed.txt";
	FILE *fp = fopen(outputFile.c_str(),"w");
	if(!fp) {
		cerr << "Error: cannot open file " << outputFile << " to save the processed data!" << endl;
		exit(-1);
	}
	
	fprintf(fp, "#Chr\tStartPos\tEndPos\tB-allele-depth\tTotal-depth\tRead-counts\tGC-content\tMappability\n");
	double *p = data_all.getEntrance();
	for(j = 0; j < T; j++) {
		for(i = 0; i < 5; i++) {
			fprintf(fp, "%d\t", (int) p[i*T+j]);
		}
		fprintf(fp, "%.2f\t", p[5*T+j]);
		fprintf(fp, "%.2f\n", p[6*T+j]);
	}
	fclose(fp);
	
	data_all.clear();
}

double GenomeData::calculateACN() {
	int i, j, chr, num_seg;
	int sindx, eindx;
	double cn, w_len = 0, t_len = 0;
	for(i = 0; i < chromosomes.size(); i++) {
		chr = chromosomes[i];
		double *p1 = data_sep[chr].getEntrance();
		Matrix<double> &segs = state_segments[chr];
		num_seg = segs.getROWS();
		for(j = 0; j < num_seg; j++) {
			sindx = segs.get(j,1);
			eindx = segs.get(j,2);
			cn = segs.get(j,3);
			w_len += cn*(p1[eindx]-p1[sindx]+1);
			t_len += p1[eindx]-p1[sindx]+1;
		}
	}
	return w_len/t_len;
}

Matrix<double>& GenomeData::getChrData(int chr) {
	if(stepsize_ds > 1) {
		return data_sep_ds[chr];
	}
	return data_sep[chr];
}

void GenomeData::setSegments(int chr, Matrix<double> &segments) {
	pthread_mutex_lock(&pm);
	state_segments[chr] = segments;
	pthread_mutex_unlock(&pm);
}

void GenomeData::fetchInputs() {
	long long spos, epos;
	string bases, qualities, sequence, ref;
	vector<long long> positions;

	//open genome reference file
	string refFile = config.getStringPara("ref");
	FastaReference fr;
    fr.open(refFile, false);
	FastaIndexEntry fie;

	//open mappability file
	string mapFile = config.getStringPara("map");
	struct bbiFile *bwf = bigWigFileOpen((char *)mapFile.c_str());
	struct bbiChromInfo *info = bbiChromList(bwf);
	struct bigWigValsOnChrom *chromVals;
	double *chromData;
	string chr_prefix;
	string name = info->name;	
	if(name.find("chrom") != string::npos) {
		chr_prefix = "chrom";
	}
	else if(name.find("chr") != string::npos) {
		chr_prefix = "chr";
	}
	else {
		chr_prefix = "";
	}
	bbiChromInfoFreeList(&info);

	//load SNP positions
	map<string, vector<long long> > positions_all;
	vector<string> chroms;
	int num_snp = readSNPs(positions_all, chroms);
	cerr << "\ntotal " << num_snp << " SNP positions are loaded from file " << config.getStringPara("snp") << "!\n" << endl;
	
	//open BAM file
	string bamFile = config.getStringPara("bam");
	BamReader reader;
	if(!reader.Open(bamFile)) {
	    cerr << "cannot open BAM file " << bamFile << endl;
	    exit(-1);
	}
	if(!reader.LocateIndex()) {
		cerr << "Building BAM index for file " << bamFile << "..." << endl;
		if(!reader.CreateIndex()) {
			cerr << "Index creation failed!" << endl;
			exit(-1);
		}
	}
	vector<RefData> refData = reader.GetReferenceData();
	string chr_prefix_bam;
	string refName = refData[0].RefName;
	if(refName.find("chrom") != string::npos) {
		chr_prefix_bam = "chrom";
	}
	else if(refName.find("chr") != string::npos) {
		chr_prefix_bam = "chr";
	}
	else {
		chr_prefix_bam = "";
	}
	
	//open output file
	string outputFile = config.getStringPara("output");
    ofstream ofs;
	ofs.open(outputFile.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << outputFile << endl;
		exit(-1);
	}
	int window = config.getIntPara("window");
	ofs << "#window size = " << window << endl;
	ofs << "Chr\tPosition\tBdepth\tTdepth\tReadCounts\tGC-content\tMappability" << endl;
	
	//struct to save data
	typedef struct snp_info{
		long long position;
		int rc, bd, td;
		char ref;
		struct snp_info *next;
	}SNPInfo;
	SNPInfo *p, *q;
	SNPInfo *rear,*head = (SNPInfo*) malloc(sizeof(SNPInfo));
	
	int is_male = config.getIntPara("male");
	int mapQ_th = config.getIntPara("mapQ_th");
	int baseQ_th = config.getIntPara("baseQ_th");
	int minDepth = config.getIntPara("minDepth");
	
	double gc, mapScore;
	long length;
	long long read_spos, read_epos;
    BamAlignment read;
    int i, j, k;
	for(i = 0; i < chroms.size(); i++) {
		string chr = chroms[i];
	    if(chr.compare("M") == 0 || chr.compare("MT") == 0) {
            continue;
        }
	    if(!is_male && (chr.compare("Y") == 0 || chr.compare("y") == 0)) {
            continue;
        }
		
		if(chr.compare("X") == 0 || chr.compare("x") == 0) {
			ref = "23";
		}
		else if(chr.compare("Y") == 0 || chr.compare("y") == 0) {
			ref = "24";
		}
		else {
			ref = chr;
		}
        
		refName = chr_prefix_bam+chr;
		for(j = 0; j < refData.size(); j++) {
			if(refName.compare(refData[j].RefName) == 0) {
				break;
			}
		}
		if(j == refData.size()) {
			//cerr << "Warning: Unrecognized chromosome name: " << chr << " from file " << snpFile << endl;
            continue;
			//exit(-1);
		}
		long refLength = refData[j].RefLength;
		int refID = reader.GetReferenceID(refName);
		
        positions = positions_all[chr];

		fie = fr.index->entry(chr);
		
        name = chr_prefix+chr;
		chromVals = bigWigValsOnChromNew();
		if(bigWigValsOnChromFetchData(chromVals, (char *) name.c_str(), bwf)){
			chromData = chromVals->valBuf;
		}
		else{
			cerr << "Can not fetch mappability data on chromosome " << chr << endl;
			exit(-1);
		}
		
		cerr << "processing reference " << chr << "...\n";
		
		head->next = NULL;
		rear = head;
		j = 0;
		k = 0;
        if(reader.SetRegion(refID, 0, refID, refLength)) {
            while(reader.GetNextAlignment(read)) {
				read_spos = read.Position+1;
				//read_epos = read.GetEndPosition(false, true)+1;
				bases = read.AlignedBases;
				read_epos = read_spos+bases.length()-1;
				qualities = read.Qualities;
                bool strand = read.IsReverseStrand();
                //string query = read.QueryBases;
                
				if(j == positions.size() && read_spos > positions[j-1]+window/2 && positions[j-1] < read_spos) {
					break;
				}
				
				if(read.IsDuplicate() || read.MapQuality < mapQ_th) {
				//if(read.MapQuality < map_quality) {
					continue;
				}

                if(bases.length() != qualities.length()) {
		    //cerr << bases << endl;
                    //cerr << qualities << endl << endl;
                    continue;
                }

				//EnQueue
				while(j < positions.size()) {
					epos = positions[j]+window/2;
					spos = epos-window+1;
					sequence = fr.getSubSequence(chr, positions[j]-1, 1);
					if(spos <= read_spos || positions[j] <= read_epos) {
						p = (SNPInfo *) malloc(sizeof(SNPInfo));
						p->position = positions[j];
						p->rc = 0;
						p->bd = 0;
						p->td = 0;
						p->ref = sequence[0];
                        if(p->ref >= 'a' && p->ref <= 'z') {
                            p->ref -= 32;
                        }
						p->next = NULL;
						rear->next = p;
						rear = p;
					}
					else {
						break;
					}
					j++;
				}
				
				p = head->next;
				while(p) {
					epos = p->position+window/2;
					spos = epos-window+1;
					//allele frequency
					if(read_spos <= p->position && read_epos >= p->position) {
						int relative_pos = p->position-read_spos;
						int q = int(qualities[relative_pos]);
						if(bases[relative_pos] != '+' && bases[relative_pos] != '-' && q >= baseQ_th) {
							if(bases[relative_pos] != p->ref) {
								p->bd++;
							}
							p->td++;
						}
					}	
					//read count
					if(read_spos >= spos && read_spos <= epos) {
						p->rc++;
					}
					if(read_spos > epos && p->position < read_spos) {	//DeQueue
						if(p->td < minDepth) {
							q = p;
							p = p->next;
							free(q);
							head->next = p;
							continue;
						}
						if(epos > fie.length) {
							epos = fie.length;
						}
						spos = epos-window+1;
						if(spos < 1) {
							spos = 1;
							epos = window;
						}
						sequence = fr.getSubSequence(chr, spos-1, window);
						calculateGC(sequence, gc, length);
						calculateMapScore(chromData, spos-1, epos-1, mapScore);
						if(gc >= 0 && mapScore <= 0.99) {
							ofs << ref << '\t' << p->position << '\t' << p->bd << '\t' << p->td << '\t' << p->rc << '\t' << gc << '\t' << mapScore << endl;
						}
						q = p;
						p = p->next;
						free(q);
						head->next = p;
					}
					else {
						p = p->next;
					}
				}
				if(head->next == NULL) {
					rear = head;
				}
				
			}
		}
		p = head->next;
		while(p) {
			if(p->td < minDepth) {
				q = p;
				p = p->next;
				free(q);
				continue;
			}
			epos = p->position+window/2;
			if(epos > fie.length) {
				epos = fie.length;
			}
			spos = epos-window+1;
			if(spos < 1) {
				spos = 1;
				epos = window;
			}
			sequence = fr.getSubSequence(chr, spos-1, window);
			calculateGC(sequence, gc, length);
			calculateMapScore(chromData, spos-1, epos-1, mapScore);
			if(gc >= 0 && mapScore <= 0.99) {
				ofs << ref << '\t' << p->position << '\t' << p->bd << '\t' << p->td << '\t' << p->rc << '\t' << gc << '\t' << mapScore << endl;
			}
			q = p;
			p = p->next;
			free(q);
		}
        bigWigValsOnChromFree(&chromVals);
	}
	
    free(head);
	ofs.close();
}

int GenomeData::readSNPs(map<string, vector<long long> > &positions_all, vector<string> &chroms) {
	int i;
	string snpFile = config.getStringPara("snp");
	FILE * fp;
	if(!(fp = fopen(snpFile.c_str(),"r"))){
		cerr << "can not open SNP file " << snpFile << endl;
		exit(-1);
	}
	int num_snp = 0;
	char buf[1000];
	char * elems[5];
	int elemnum;
	long line_num = 0;
    while(fgets(buf, 1000, fp)){
		line_num++;
		if(buf[0] == '#') {
			continue;
		}
		elemnum = split(buf, '\t', elems);
		if(elemnum != 2) {
			cerr << "Error: malformed SNP file " << snpFile <<
                    ", there should be 2 fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		num_snp++;
		string chr = elems[0]; 
		chr = abbrOfChr(chr);
		long long pos = atoll(elems[1]);
		positions_all[chr].push_back(pos);
		for(i = 0; i < chroms.size(); i++) {
			if(chr.compare(chroms[i]) == 0) {
				break;
			}
		}
		if(i == chroms.size()) {
			chroms.push_back(chr);
		}
	}
	
	fclose(fp);
	return num_snp;
}

void GenomeData::calculateGC(string sequence, double &gc_content, long &length) {
	long gc_count = 0;
	long n_count = 0;
	char c;
	length = 0;

	for(int i = 0; i < sequence.length(); i++){
		c = sequence.at(i);
		if(c == '>'){
			break;
		}
		switch(c) {
			case 'G':
			case 'C':
			case 'g':
			case 'c':
				gc_count++;
				length++;
				break;
			case 'A':
			case 'T':
			case 'a':
			case 't':
				length++;
				break;
			default:
				n_count++;
				length++;
				break;
		}
	}

	gc_content = (n_count < length)? (double) gc_count/(length-n_count) : -1;
}

void GenomeData::calculateMapScore(double * chromData, long long spos, long long epos, double &mapScore) {
	mapScore = 0;
	for(int i = spos; i <= epos; i++){
		mapScore += chromData[i];
	}
	mapScore = mapScore/(epos-spos+1);
}


