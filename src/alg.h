#pragma once
#include <string>
#include "Preprocess.h"
#include "mf_alsh.h"
#include "performance.h"
#include "basis.h"
#include "method.h"
#include "linscan.h"
#include "wand.h"
#include "sinnamon.h"
#include "sosia.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <cstring>
#include <chrono>

struct resOutput
{
	std::string algName;
	int L;
	int K;
	int k;
	float ub;
	float time;
	float recall;
	float ratio;
	float cost;
	float kRatio;
};

extern std::string data_fold, index_fold;
extern std::string data_fold1, data_fold2;

#if defined(unix) || defined(__unix__)
inline void localtime_s(tm* ltm, time_t* now) {}
#endif

inline resOutput FargoSearch(fargo::Hash& myslsh, float c_, int m_, int k_, int L_, int K_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	lsh::progress_display pd(Qnum);
	Performance<fargo::Query> perform;
	lsh::timer timer1;
#pragma omp parallel for
	for (int j = 0; j < Qnum; j++){
		fargo::Query query(j, c_, k_, myslsh, prep, m_);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "FARGO";
	res.L = myslsh.L;
	res.K = myslsh.K;
	res.ub = m_;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)myslsh.N);
	res.kRatio=perform.kRatio/perform.num;
	//delete[] ltm;
	return res;
}

inline resOutput linscanSearch(Linscan& lins, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	if (Qnum > prep.queries->n)Qnum = prep.queries->n;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
//#pragma omp parallel for num_threads(32)
	for (int j = 0; j < Qnum; j++){
		Query query(j, c_, k_, prep, ub_);
		lins.knnSearch(query);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "PRODUCT COST:      " << perform.cnt_product << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "Linscan";
	res.k = k_;
	res.L = -1;
	res.K = -1;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)lins.n);
	res.kRatio = perform.kRatio / perform.num;
	//delete[] ltm;
	return res;
}

inline resOutput wandSearch0(wand& wand, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	if (Qnum > prep.queries->n)Qnum = prep.queries->n;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
	//#pragma omp parallel for num_threads(32)
	for (int j = 0; j < Qnum; j++) {
		Query query(j, c_, k_, prep, ub_);
		wand.knnSearch0(query);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "PRODUCT COST:      " << perform.cnt_product << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "wand0";
	res.k = k_;
	res.L = -1;
	res.K = -1;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)wand.n);
	res.kRatio = perform.kRatio / perform.num;
	//delete[] ltm;
	return res;
}

inline resOutput wandSearch1(wand& wand, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	if (Qnum > prep.queries->n)Qnum = prep.queries->n;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
	//#pragma omp parallel for num_threads(32)
	for (int j = 0; j < Qnum; j++) {
		Query query(j, c_, k_, prep, ub_);
		wand.knnSearch1(query);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "PRODUCT COST:      " << perform.cnt_product << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "wand1";
	res.k = k_;
	res.L = -1;
	res.K = -1;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)wand.n);
	res.kRatio = perform.kRatio / perform.num;
	//delete[] ltm;
	return res;
}

inline resOutput wandSearch2(wand& wand, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	if (Qnum > prep.queries->n)Qnum = prep.queries->n;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
	//#pragma omp parallel for num_threads(32)
	for (int j = 0; j < Qnum; j++) {
		Query query(j, c_, k_, prep, ub_);
		wand.knnSearch2(query);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "PRODUCT COST:      " << perform.cnt_product << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "wand2";
	res.L = -1;
	res.K = -1;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)wand.n);
	res.kRatio = perform.kRatio / perform.num;
	//delete[] ltm;
	return res;
}

inline resOutput wandSearch3(wand& wand, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	if (Qnum > prep.queries->n)Qnum = prep.queries->n;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
	//#pragma omp parallel for num_threads(32)
	for (int j = 0; j < Qnum; j++) {
		Query query(j, c_, k_, prep, ub_);
		wand.knnSearch3(query);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "PRODUCT COST:      " << perform.cnt_product << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "wand3";
	res.L = -1;
	res.K = -1;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)wand.n);
	res.kRatio = perform.kRatio / perform.num;
	//delete[] ltm;
	return res;
}

inline resOutput sinaSearch(Sinnamon& sina, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	if (Qnum > prep.queries->n)Qnum = prep.queries->n;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
	//#pragma omp parallel for num_threads(32)
	for (int j = 0; j < Qnum; j++) {
		Query query(j, c_, k_, prep, ub_);
		sina.knnSearch(query);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "PRODUCT COST:      " << perform.cnt_product << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "Sinnamon";
	res.k = k_;
	res.L = -1;
	res.K = -1;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)sina.n);
	res.kRatio = perform.kRatio / perform.num;
	//delete[] ltm;
	return res;
}

inline resOutput bruteforce(Linscan& lins, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
//#pragma omp parallel for
	for (int j = 0; j < Qnum; j++){
		Query query(j, c_, k_, prep, ub_);
		lins.bruteforce(query, prep.data);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "bruteforce";
	res.L = -1;
	res.K = -1;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)lins.n);
	res.kRatio = perform.kRatio / perform.num;
	return res;
}

inline resOutput ivf_BF(ivfForSet& sets, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
	for (int j = 0; j < Qnum; j++)
	{
		Query query(j, c_, k_, prep, ub_);
		sets.knnSearch1(query);
		perform.update(query, prep);
		++pd;
	}

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = "ivf_BF";
	res.k = k_;
	res.L = -1;
	res.K = mips2set::l;
	res.ub = -1;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)sets.n);
	res.kRatio = perform.kRatio / perform.num;
	return res;
}

template <class myMinHash>
inline resOutput minHash(myMinHash& mh, float c_, int k_, int ub_, Preprocess& prep)
{
	std::string query_result = ("results/MF_ALSH_result.csv");

	lsh::timer timer;
	std::cout << std::endl << "RUNNING QUERY ..." << std::endl;

	int Qnum = 100;
	lsh::progress_display pd(Qnum);
	Performance<Query> perform;
	lsh::timer timer1;
	for(int t=0;t<1;t++){
		for (int j = 0; j < Qnum; j++){
			Query query(j, c_, k_, prep, ub_);
			mh.knnSearch2(query);
			perform.update(query, prep);
			++pd;
		}
	}
	

	float mean_time = (float)perform.time_total / perform.num;
	std::cout << "AVG HASH TIME:     " << (float)perform.time_hash / perform.num * 1000 << "ms." << std::endl;
	std::cout << "AVG SIFT TIME:     " << (float)perform.time_sift / perform.num * 1000 << "ms." << std::endl;
	std::cout << "AVG VERIFY TIME:   " << (float)perform.time_verify / perform.num * 1000 << "ms." << std::endl;
	std::cout << "AVG QUERY TIME:    " << mean_time * 1000 << "ms." << std::endl;
	std::cout << "PRODUCT COST:      " << perform.cnt_product << std::endl;
	std::cout << "AVG CNT RATIO:    " << (float)perform.cnt / perform.num / prep.data->n << std::endl;
	std::cout << "AVG RECALL:        " << ((float)perform.NN_num) / (perform.num * k_) << std::endl;
	std::cout << "AVG RATIO:         " << ((float)perform.ratio) / (perform.res_num) << std::endl << std::endl;

	time_t now = time(0);
	tm* ltm = new tm[1];
	localtime_s(ltm, &now);

	resOutput res;
	res.algName = mh.getAlgName();
	res.k = k_;
	res.L = mh.r;
	res.K = mips2set::l;// static member
	res.ub = ub_;
	res.time = mean_time * 1000;
	res.recall = ((float)perform.NN_num) / (perform.num * k_);
	res.ratio = ((float)perform.ratio) / (perform.res_num);
	res.cost = ((float)perform.cost) / ((long long)perform.num);
	//res.cost = ((float)perform.cost) / ((long long)perform.num * (long long)prep.data->n);
	res.kRatio = perform.kRatio / perform.num;
	return res;
}

void saveAndShow(float c, int k, std::string& dataset, std::vector<resOutput>& res);


