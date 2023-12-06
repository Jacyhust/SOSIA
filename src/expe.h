#pragma once
#include "alg.h"

void runAll(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {
	Linscan lins(prep.data);
	res.push_back(linscanSearch(lins, c, k, ub, prep));
	// res.push_back(bruteforce(lins, c, k, ub, prep));
	// return 0;

	wand wand(lins, prep.data, F);
	// res.push_back(wandSearch0(wand, c, k, 0, prep));
	res.push_back(wandSearch1(wand, c, k, 0, prep));
	// res.push_back(wandSearch2(wand, c, k, 0, prep));
	// res.push_back(wandSearch3(wand, c, k, 0, prep));

	//return 0;

	mips2set sets(prep, l);
	ivfForSet ivf(sets, prep);
	res.push_back(ivf_BF(ivf, c, k, ub, prep));

	//myMinHashv2 mhash(prep, sets, b, r);
	//res.push_back(minHash(mhash, c, k, ub, prep));

	myMinHashBase minBase(prep, sets, 1, m);

	//Sosia sosia(prep, sets, b, r);
	Sosia sosia(minBase);
	res.push_back(minHash(sosia, c, k, ub * 2, prep));

	//myMinHashv3 mhashv3(prep, sets, b, r);
	// myMinHashv3 mhashv3(minBase);
	// auto minRes = minHash(mhashv3, c, k, ub * 2, prep);
	// res.push_back(minRes);

	auto res_sosia=res[res.size()-1];
	float T = 2 * res_sosia.time / 1000;
	std::cout << "T=" << T << std::endl;
	Sinnamon sina(lins, prep.data, sin_m, sin_h, T);
	res.push_back(sinaSearch(sina, c, k, ub / 2, prep));
}

void varyk(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {

	std::vector<int> ks={1,10,20,30,40,50,60,70,80,90,100};

	Linscan lins(prep.data);
	//for(auto& k:ks) res.push_back(linscanSearch(lins, c, k, ub, prep));
	// res.push_back(bruteforce(lins, c, k, ub, prep));
	// return 0;

	wand wand(lins, prep.data, F);
	// res.push_back(wandSearch0(wand, c, k, 0, prep));
	//for(auto& k:ks) res.push_back(wandSearch1(wand, c, k, 0, prep));
	// res.push_back(wandSearch2(wand, c, k, 0, prep));
	// res.push_back(wandSearch3(wand, c, k, 0, prep));

	//return 0;

	mips2set sets(prep, l);
	ivfForSet ivf(sets, prep);
	//for(auto& k:ks) res.push_back(ivf_BF(ivf, c, k, ub, prep));

	//myMinHashv2 mhash(prep, sets, b, r);
	//res.push_back(minHash(mhash, c, k, ub, prep));

	myMinHashBase minBase(prep, sets, 1, m);

	//Sosia sosia(prep, sets, b, r);
	Sosia sosia(minBase);
	//for(auto& k:ks) res.push_back(minHash(sosia, c, k, ub * 2, prep));
	res.push_back(minHash(sosia, c, k, ub * 2, prep));
	//myMinHashv3 mhashv3(prep, sets, b, r);
	// myMinHashv3 mhashv3(minBase);
	// auto minRes = minHash(mhashv3, c, k, ub * 2, prep);
	// res.push_back(minRes);

	auto res_sosia=res[res.size()-1];
	float T = 2 * res_sosia.time / 1000;
	std::cout << "T=" << T << std::endl;
	Sinnamon sina(lins, prep.data, sin_m, sin_h, T);
	for(auto& k:ks) res.push_back(sinaSearch(sina, c, k, ub / 2, prep));
}

void varyL(Preprocess& prep, int c, int k, int r, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {
	std::vector<int> ls = { 1,10,20,40,80 };
	for (auto& l : ls) {
		mips2set sets(prep, l);
		// ivfForSet ivf(sets, prep);
		// res.push_back(ivf_BF(ivf, c, k, ub, prep));
		myMinHashBase minBase(prep, sets, 1, r);

		//Sosia sosia(prep, sets, b, r);
		Sosia sosia(minBase);
		res.push_back(minHash(sosia, c, k, ub * 2, prep));

		//myMinHashv3 mhashv3(prep, sets, b, r);
		// myMinHashv3 mhashv3(minBase);
		// auto minRes = minHash(mhashv3, c, k, ub * 2, prep);
		// res.push_back(minRes);
	}
}


void varyM(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {
	//std::vector<int> ms = { 25,50,100,150,200 };
	std::vector<int> ms = { 25,50,100,150,200};
	mips2set sets(prep, l);
	// ivfForSet ivf(sets, prep);
	// res.push_back(ivf_BF(ivf, c, k, ub, prep));
	for (auto& m : ms) {
		myMinHashBase minBase(prep, sets, 1, m);

		Sosia sosia(minBase);
		res.push_back(minHash(sosia, c, k, ub * 2, prep));

		//myMinHashv3 mhashv3(prep, sets, b, r);
		// myMinHashv3 mhashv3(minBase);
		// auto minRes = minHash(mhashv3, c, k, ub * 2, prep);
		// res.push_back(minRes);
	}
}

void RecallTime(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {

	//std::vector<int> ks={1,10,20,30,40,50,60,70,80,90,100};

	Linscan lins(prep.data);
	
	res.push_back(linscanSearch(lins, c, k, ub, prep));
	// res.push_back(bruteforce(lins, c, k, ub, prep));
	// return 0;
	// std::vector<float> Fs={1.0,2.0,4.0,8.0,16.0,32.0};
	std::vector<float> Fs={1.1,1.2,1.4,1.6,1.8};
	wand wand(lins, prep.data, F);
	// res.push_back(wandSearch0(wand, c, k, 0, prep));
	for(auto& F:Fs){
		wand.F=F;
		res.push_back(wandSearch1(wand, c, k, 0, prep));
	} 
	// res.push_back(wandSearch2(wand, c, k, 0, prep));
	// res.push_back(wandSearch3(wand, c, k, 0, prep));

	//return 0;

	mips2set sets(prep, l);
	ivfForSet ivf(sets, prep);
	res.push_back(ivf_BF(ivf, c, k, ub, prep));

	myMinHashBase minBase(prep, sets, 1, m);
	Sosia sosia(minBase);
	//std::vector<int> qms={20,60,100,140,180,220,280,340,400};
	// std::vector<int> ts={5000,10000,15000,20000,40000,80000};
	// for(auto& ub:ts) res.push_back(minHash(sosia, c, k, ub * 2, prep));
	res.push_back(minHash(sosia, c, k, ub * 2, prep));
	// //myMinHashv3 mhashv3(prep, sets, b, r);
	// // myMinHashv3 mhashv3(minBase);
	// // auto minRes = minHash(mhashv3, c, k, ub * 2, prep);
	// // res.push_back(minRes);

	//varyM(prep, c, k, m, ub, F, l, sin_m, sin_h, res);

	auto res_sosia=res[res.size()-1];
	float T = 2 * res_sosia.time / 1000;
	std::cout << "T=" << T << std::endl;
	Sinnamon sina(lins, prep.data, sin_m, sin_h, T);

	

	std::vector<float> Ts={0.25,0.5,1.0,2.0,4.0,8.0,16.0,32.0};
	for(auto& t:Ts){
		sina.T=0.5*t*T;
		res.push_back(sinaSearch(sina, c, k, ub / 2, prep));
	} 
	
}