#pragma once
#include <string>
#include "Preprocess.h"
#include "performance.h"
#include "basis.h"
#include "method.h"
//#include "linscan.h"
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


