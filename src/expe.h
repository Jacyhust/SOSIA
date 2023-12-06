#pragma once
#include "alg.h"

void runAll(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {
	mips2set sets(prep, l);
	ivfForSet ivf(sets, prep);
	myMinHashBase minBase(prep, sets, 1, m);
	Sosia sosia(minBase);
	res.push_back(minHash(sosia, c, k, ub * 2, prep));
}

void varyk(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {

	std::vector<int> ks={1,10,20,30,40,50,60,70,80,90,100};

	mips2set sets(prep, l);
	ivfForSet ivf(sets, prep);
	myMinHashBase minBase(prep, sets, 1, m);
	Sosia sosia(minBase);
	for(auto& k:ks) res.push_back(minHash(sosia, c, k, ub * 2, prep));
}

void varyL(Preprocess& prep, int c, int k, int r, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {
	std::vector<int> ls = { 1,10,20,40,80 };
	for (auto& l : ls) {
		mips2set sets(prep, l);
		myMinHashBase minBase(prep, sets, 1, r);
		Sosia sosia(minBase);
		res.push_back(minHash(sosia, c, k, ub * 2, prep));
	}
}


void varyM(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {
	std::vector<int> ms = { 25,50,100,150,200};
	mips2set sets(prep, l);
	for (auto& m : ms) {
		myMinHashBase minBase(prep, sets, 1, m);
		Sosia sosia(minBase);
		res.push_back(minHash(sosia, c, k, ub * 2, prep));
	}
}

void RecallTime(Preprocess& prep, int c, int k, int m, int ub, int F, int l, int sin_m, int sin_h, std::vector<resOutput>& res) {

	mips2set sets(prep, l);
	ivfForSet ivf(sets, prep);
	myMinHashBase minBase(prep, sets, 1, m);
	Sosia sosia(minBase);
	std::vector<int> ts={5000,10000,15000,20000,40000,80000};
	for(auto& ub:ts) res.push_back(minHash(sosia, c, k, ub * 2, prep));
}