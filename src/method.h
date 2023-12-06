#pragma once
#include "def.h"
#include "Preprocess.h"
#include "basis.h"
#include <cmath>
#include <assert.h>
#include <vector>
#include <queue>
#include <cfloat>
#include <random>
#include <ctime>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
using trans_set = std::vector<int>*;

class mips2set
{
public:
	static int l; //default length per dim.
	size_t n = 0;
	int dim = 0;
	float maxv = -FLT_MAX;
	static std::mt19937 rng;
	static std::uniform_real_distribution<float> ur;

	trans_set transSets = nullptr;
public:
	mips2set(Preprocess& prep, int l_=10) {
		mips2set::l = l_;
		n = prep.data->n;
		dim = prep.data->dim;
		maxValue(prep.data);
		std::cout << "SOS TRANSFORMING ...\n";
		lsh::timer timer;
		getSets(prep.data);
		std::cout << "TRANSFORMING TIME: " << timer.elapsed() << " s.\n\n";
		
	}

	void maxValue(SparseData* data) {
		for (size_t i = 0; i < data->nnz; ++i) {
			if (maxv < data->val[i]) {
				maxv = data->val[i];
			}
		}
	}

	void getSets(SparseData* data) {
		transSets = new std::vector<int>[n];

		lsh::progress_display pd(n);
#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			auto point = data->generateSV(i);
			for (int j = 0; j < point->nnz; ++j) {
				float pthred = point->val[j] / maxv;
				for (int r = point->indices[j] * l; r < point->indices[j] * l + l; ++r) {
					if (ur(rng) < pthred) transSets[i].push_back(r);
				}
			}
			++pd;
		}
	}

	static void getQuerySet(Query& q, std::vector<int>& qset) {
		auto point = &(q.queryPoint);
		float qmaxv = point->getMax();
		for (int j = 0; j < point->nnz; ++j) {
			float pthred = point->val[j] / qmaxv;
			for (int r = point->indices[j] * l; r < point->indices[j] * l + l; ++r) {
				if (ur(rng) < pthred) qset.push_back(r);
			}
		}
	}

	// static void getQuerySet(Query& q, std::vector<int>& qset) {
	// 	auto point = &(q.queryPoint);
	// 	float qmaxv = point->getMax();
	// 	for (int j = 0; j < point->nnz; ++j) {
	// 		float pthred = point->val[j] / qmaxv;
	// 		pthred*=l;
	// 		int num=round(pthred);
	// 		for (int r = point->indices[j] * l; r < point->indices[j] * l + num; ++r) {
	// 			qset.push_back(r);
	// 		}
	// 	}
	// }
};

class ivfForSet
{
	std::vector<std::vector<size_t>> ivfIndex;
public:
	size_t n = 0;
	int dim = 0;
	int maxElem = -1;
	ivfForSet() = default;
	SparseData* data = nullptr;
	ivfForSet(mips2set& sets, Preprocess& prep) {
		n = sets.n;
		dim = sets.dim;
		maxElem = dim * sets.l;
		data = prep.data;

		std::cout << "IVF4SET INDEXING ...\n";
		lsh::timer timer;
		buildIndex(sets);
		std::cout << "INDEXING TIME: " << timer.elapsed() << " s.\n";
	}

	void buildIndex(mips2set& sets) {
		ivfIndex.resize(maxElem);
		lsh::progress_display pd(n);
		for (size_t i = 0; i < n; ++i) {
			for (auto& id : sets.transSets[i])
				ivfIndex[id].emplace_back(i);

			++pd;
		}
	}

	void knnSearch(Query& q) {
		lsh::timer timer;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset);

		std::vector<int> score(n, 0);
		for (auto& ind : qset) {
			for (auto& item : ivfIndex[ind]) {
				score[item] ++;
			}
		}

		q.time_total = timer.elapsed();
		return;

		int ub = 200;
		//ub = q.k;
		Res* res_PQ = new Res[ub + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			res_PQ[size] = Res(i, score[i]);
			if (size < ub) {
				size++;
				std::push_heap(res_PQ, res_PQ + size);
			}
			else if (res_PQ[0].inp < res_PQ[size].inp) {
				size++;
				std::push_heap(res_PQ, res_PQ + size);
				std::pop_heap(res_PQ, res_PQ + size);
				size--;
			}
		}
		  
		//Linscan lins(data, res_PQ, size);
		//lins.knnSearch(q);
		//for (int i = 0; i < q.k; ++i) q.res[i].id = res_PQ[q.res[i].id].id;

		for (int i = 0; i < size; ++i) {
			res_PQ[i].inp = q.queryPoint.dotProduct2SV(data->generateSV(res_PQ[i].id));
			//res_PQ[i].inp = data->dotProduct2SV(res_PQ[i].id, q.qvec);
		}
		Res* res1 = new Res[q.k + 1];
		int size1 = 0;
		for (int i = 0; i < size; ++i) {
			res1[size1] = res_PQ[i];
			if (size1 < q.k) {
				size1++;
				std::push_heap(res1, res1 + size1);
			}
			else if (res1[0].inp < res1[size1].inp) {
				size1++;
				std::push_heap(res1, res1 + size1);
				std::pop_heap(res1, res1 + size1);
				size1--;
			}
		}

		int len = size1;
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr] = res1[0];
			std::pop_heap(res1, res1 + size1);
			size1--;
			rr--;
		}

		std::vector<Res> res(q.k + 1);
		q.time_total = timer.elapsed();
	}

	void knnSearch1(Query& q) {
		lsh::timer timer;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset);

		std::vector<ResInt> score(n);
		for (int i = 0; i < n; ++i) score[i].id = i;

		for (auto& ind : qset) {
			for (auto& item : ivfIndex[ind]) {
				score[item].inp++;
			}
		}
		//int ub = 200;
		int ub = q.k + q.ub;
		if (ub > n) ub = n;
		std::partial_sort(score.begin(), score.begin() + ub, score.end());
		//for (int j = 0; j < q.queryPoint.nnz; ++j) {
		//	int ind = q.queryPoint.indices[j];
		//	float val = q.queryPoint.val[j];
		//	for (auto& item : ivfIndex[ind]) {
		//		score[item.id] += val * item.val;
		//	}
		//}
		//q.time_total = timer.elapsed();
		//if (ub > q.res.size()) ub = q.res.size();
		q.res.resize(ub);
		for (int i = 0; i < ub; ++i) q.res[i] = Res(score[i].id, data->dotProduct2SV(score[i].id, q.qvec));
		std::partial_sort(q.res.begin(), q.res.begin() + q.k, q.res.end());
		q.res.resize(q.k);
		q.time_total = timer.elapsed();
		return;

		
		Res* res_PQ = new Res[ub + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			res_PQ[size] = Res(i, score[i].inp);
			if (size < ub) {
				size++;
				std::push_heap(res_PQ, res_PQ + size);
			}
			else if (res_PQ[0].inp < res_PQ[size].inp) {
				size++;
				std::push_heap(res_PQ, res_PQ + size);
				std::pop_heap(res_PQ, res_PQ + size);
				size--;
			}
		}

		//Linscan lins(data, res_PQ, size);
		//lins.knnSearch(q);
		//for (int i = 0; i < q.k; ++i) q.res[i].id = res_PQ[q.res[i].id].id;

		for (int i = 0; i < size; ++i) {
			res_PQ[i].inp = q.queryPoint.dotProduct2SV(data->generateSV(res_PQ[i].id));
			//res_PQ[i].inp = data->dotProduct2SV(res_PQ[i].id, q.qvec);
		}
		Res* res1 = new Res[q.k + 1];
		int size1 = 0;
		for (int i = 0; i < size; ++i) {
			res1[size1] = res_PQ[i];
			if (size1 < q.k) {
				size1++;
				std::push_heap(res1, res1 + size1);
			}
			else if (res1[0].inp < res1[size1].inp) {
				size1++;
				std::push_heap(res1, res1 + size1);
				std::pop_heap(res1, res1 + size1);
				size1--;
			}
		}

		int len = size1;
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr] = res1[0];
			std::pop_heap(res1, res1 + size1);
			size1--;
			rr--;
		}

		q.time_total = timer.elapsed();
	}
};

class myMinHash
{
	int n = 0;
	int b = 0;
	int r = 0;
	int s = 0;
	int len = 0;
	mips2set* mysets = nullptr;
	SparseData* data = nullptr;
	int** hashFuns = nullptr;
	int** hashvals = nullptr;
	int bucketSize = 0;
	std::vector<std::unordered_multimap<int, int>> myIndexes;
public:
	myMinHash(Preprocess& prep, mips2set& sets, int b_, int r_) {
		n = sets.n;
		b = b_;
		r = r_;
		s = b * r;
		len = sets.dim * sets.l;
		bucketSize = (1 << 10) - 1;
		bucketSize = len;
		data = prep.data;
		mysets = &sets;
		setHash();
		getHash(sets);
		buildIndex();
	}

	void setHash() {
		hashFuns = new int* [s];
		for (int i = 0; i < s; ++i) {
			hashFuns[i] = new int[len];
			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
			for (int j = 0; j < len; ++j) hashFuns[i][j] = j;
			std::shuffle(hashFuns[i], hashFuns[i] + len, std::default_random_engine(seed));
		}
	}

	void getHash(mips2set& sets) {
		hashvals= new int* [n];
		for (int i = 0; i < n; ++i)hashvals[i] = new int[s];
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < s; ++j) hashvals[i][j] = len;
			auto point = sets.transSets[i];
			for (auto& id : point) {
				for (int j = 0; j < s; ++j) if (hashvals[i][j] > hashFuns[j][id]) hashvals[i][j] = hashFuns[j][id];
			}
		}
	}

	int getKeys(int* hs) {
		size_t res = 0;
		for (int j = 0; j < b; ++j) {
			res *= (len % bucketSize);
			//res *= len;
			res %= bucketSize;
			//res = fastrange64(res * len, bucketSize);
			res += hs[j] % bucketSize;
			//res += fastrange64(hs[j], bucketSize);
		}
		return res;
	}

	void buildIndex() {
		//int width = floor(pow(bucketSize, 1.0 / b));
		myIndexes.resize(r);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < r; ++j) {
				int key = getKeys(hashvals[i] + j * b);
				myIndexes[j].insert({ key,i });
			}
		}


	}

	void knnSearch(Query& q) {
		lsh::timer timer;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset);
		std::vector<int> qhash(s, len);
		std::vector<bool> unseen(n, true);
		for (auto& id : qset) {
			for (int j = 0; j < s; ++j) if (qhash[j] > hashFuns[j][id]) qhash[j] = hashFuns[j][id];
		}

		int ub = 200;
		Res* cands = new Res[ub];
		int cnt = 0;
		for (int j = 0; j < r; ++j) {
			int key = getKeys(&(qhash[0]) + j * b);
			auto pr = myIndexes[j].equal_range(key);
			while (pr.first != pr.second){
				if (unseen[pr.first->second]){
					//res_pair.id = pr.first->second;
					int id = pr.first->second;
					float inp = q.queryPoint.dotProduct2SV(data->generateSV(id));
					cands[cnt++] = Res(id, inp);
					unseen[pr.first->second] = false;
					if (cnt == ub) break;
				}
				++pr.first; // Increment begin iterator
			}

			if (cnt == ub) break;
		}

		std::sort(cands, cands + cnt);
		q.cost = cnt;
		if (q.cost <= q.k)q.res.assign(cands, cands + cnt);
		else q.res.assign(cands, cands + q.k);

		q.time_total = timer.elapsed();
	}

	void knnSearch1(Query& q) {
		lsh::timer timer;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset);
		std::vector<int> qhash(s, len);
		std::vector<bool> unseen(n, true);
		for (auto& id : qset) {
			for (int j = 0; j < s; ++j) if (qhash[j] > hashFuns[j][id]) qhash[j] = hashFuns[j][id];
		}

		std::unordered_map<int, int> candRes;
		for (int j = 0; j < r; ++j) {
			int key = getKeys(&(qhash[0]) + j * b);
			auto pr = myIndexes[j].equal_range(key);
			while (pr.first != pr.second) {
				int id = pr.first->second;
				if (candRes.find(id) != candRes.end()) candRes[id]++;
				else candRes.insert({ id,1 });

				++pr.first; // Increment begin iterator
			}
		}

		int size = candRes.size();
		Res* cands = new Res[size];
		int cnt = 0;
		for (auto iter = candRes.begin(); iter != candRes.end(); iter++) {
			int id = iter->first;
			float jaccord_inv = (float)r / iter->second;
			float overlap = (float)(qset.size() + mysets->transSets[id].size());
			overlap /= 1 + jaccord_inv;
			cands[cnt++] = Res(iter->first, overlap);
			//cands[cnt++] = Res(iter->first, iter->second);
		}
		int ub = 2000;
		ub = q.k + 200;
		//  ub = size;
		if (ub > cnt) {
			q.res.assign(cands, cands + cnt);
			ub = cnt;
		}
		else {
			std::partial_sort(cands, cands + ub, cands + cnt);
			
		}
		q.res.resize(ub);
		
		for (int i = 0; i < ub; ++i) q.res[i] = Res(cands[i].id, data->dotProduct2SV(cands[i].id, q.qvec));
		std::partial_sort(q.res.begin(), q.res.begin() + q.k, q.res.end());
		q.res.resize(q.k);
		q.time_total = timer.elapsed();
		return; 
		q.time_total = timer.elapsed();
	}

	void knnSearch2(Query& q) {
		lsh::timer timer;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset);
		std::vector<int> qhash(s, len);
		std::vector<bool> unseen(n, true);
		for (auto& id : qset) {
			for (int j = 0; j < s; ++j) if (qhash[j] > hashFuns[j][id]) qhash[j] = hashFuns[j][id];
		}

		std::vector<int> candRes(n,0);
		for (int j = 0; j < r; ++j) {
			int key = getKeys(&(qhash[0]) + j * b);
			auto pr = myIndexes[j].equal_range(key);
			while (pr.first != pr.second) {
				candRes[pr.first->second]++;
				//int id = pr.first->second;
				//if (candRes.find(id) != candRes.end()) candRes[id]++;
				//else candRes.insert({ id,1 });

				++pr.first; // Increment begin iterator
			}
		}

		int ub = 2000;
		ub = q.k + 20000; 
		Res* cands = new Res[ub + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			if (candRes[i] > 0) {
				float jaccord_inv = (float)r / candRes[i];
				float overlap = (float)(qset.size() + mysets->transSets[i].size());
				overlap /= 1 + jaccord_inv;
				cands[size] = Res(i, overlap);
				//cands[size] = Res(i, candRes[i]);
				if (size < ub) {
					size++;
					std::push_heap(cands, cands + size);
				}
				else if (cands[0].inp < cands[size].inp) {
					size++;
					std::push_heap(cands, cands + size);
					std::pop_heap(cands, cands + size);
					size--;
				}
			}
			
		}
		if (size < ub) ub = size;
		q.res.resize(ub);

		for (int i = 0; i < ub; ++i) q.res[i] = Res(cands[i].id, data->dotProduct2SV(cands[i].id, q.qvec));
		std::partial_sort(q.res.begin(), q.res.begin() + q.k, q.res.end());
		q.res.resize(q.k);
		q.time_total = timer.elapsed();
		return;
	}
};

class myMinHashBase
{
public:
	int n = 0;
	int b = 0;
	int r = 0;
	int s = 0;
	int len = 0;
	mips2set* mysets = nullptr;
	SparseData* data = nullptr;
	int** hashFuns = nullptr;
	int** hashvals = nullptr;
public:
	myMinHashBase(Preprocess& prep, mips2set& sets, int b_, int r_) {
		n = sets.n;
		b = b_;
		r = r_;
		s = b * r;
		len = sets.dim * sets.l + 1;//The last bucket is set for the empty set
		data = prep.data;
		mysets = &sets;
		lsh::timer timer;
		std::cout << "SETTING HASH PARAMETER..." << std::endl;
		timer.restart();
		setHash();
		std::cout << "SETTING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

		std::cout << "COMPUTING HASH..." << std::endl;
		timer.restart();
		getHash(sets);
		std::cout << "COMPUTING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

	}

	void setHash() {
		hashFuns = new int* [s];
		lsh::progress_display pd(s);
#pragma omp parallel for
		for (int i = 0; i < s; ++i) {
			hashFuns[i] = new int[len];
			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
			for (int j = 0; j < len; ++j) hashFuns[i][j] = j;
			std::shuffle(hashFuns[i], hashFuns[i] + len, std::default_random_engine(seed));
			++pd;
		}
	}

	void getHash(mips2set& sets) {
		hashvals = new int* [n];
		for (int i = 0; i < n; ++i)hashvals[i] = new int[s];

		lsh::progress_display pd(n);
#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < s; ++j) hashvals[i][j] = len - 1;
			auto point = sets.transSets[i];
			for (auto& id : point) { // Incur an error when the set is empty
				for (int j = 0; j < s; ++j) if (hashvals[i][j] > hashFuns[j][id]) hashvals[i][j] = hashFuns[j][id];
			}
			++pd;
		}
	}

	//virtual void knnSearch2(Query& q) = 0;
};

class myMinHashv2 :public myMinHashBase
{
	const std::string algName = "minHashv2";
	int bucketSize = 0;
	std::vector<std::unordered_multimap<int, int>> myIndexes;
public:
	myMinHashv2(Preprocess& prep, mips2set& sets, int b_, int r_) :myMinHashBase(prep, sets, b_, r_) {
		bucketSize = (1 << 10) - 1;
	
		lsh::timer timer;
		std::cout << "\nMinv2 INDEXING..." << std::endl;
		timer.restart();
		buildIndex();
		std::cout << "BUILDING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;
	}

	std::string getAlgName() const {
		return algName;
	}

	int getKeys(int* hs) {
		size_t res = 0;
		for (int j = 0; j < b; ++j) {
			res *= (len % bucketSize);
			//res *= len;
			res %= bucketSize;
			//res = fastrange64(res * len, bucketSize);
			res += hs[j] % bucketSize;
			//res += fastrange64(hs[j], bucketSize);
		}
		return res;
	}

	void buildIndex() {
		myIndexes.resize(r);
		lsh::progress_display pd(r);
#pragma omp parallel for
		for (int j = 0; j < r; ++j) {
			for (int i = 0; i < n; ++i) {
				for (int k = 0; k < b; ++k) {
					int key = hashvals[i][j * b + k];
					myIndexes[j].insert({ key,i });
				}
			}
			++pd;
		}
	}

	void knnSearch2(Query& q) {
		lsh::timer timer, timer1;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset);
		std::vector<int> qhash(s, len);
		std::vector<bool> unseen(n, true);
		for (auto& id : qset) {
			for (int j = 0; j < s; ++j) if (qhash[j] > hashFuns[j][id]) qhash[j] = hashFuns[j][id];
		}
		
		std::vector<int> candRes(n, 0);
		size_t cnt = 0;
		for (int j = 0; j < r; ++j) {
			std::unordered_set<int> cands;
			for (int k = 0; k < b; ++k) {
				int key = qhash[j * b + k];
				auto pr = myIndexes[j].equal_range(key);
				while (pr.first != pr.second) {
					candRes[pr.first->second]++;
					//cands.insert(pr.first->second);
					++pr.first; // Increment begin iterator
					cnt++;
				}
			}
		}
		q.cnt = cnt;
		q.time_hash = timer1.elapsed();
		timer1.restart();

		//int ub = 2e3;
		int ub = q.k + q.ub;
		Res* cands = new Res[ub + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			if (candRes[i] > 0) {
				float jaccord_inv = (float)r / candRes[i];
				//float jaccord_inv = 1.0 / (1 - std::pow(1 - candRes[i] / (float)r, 1.0 / b));
				float overlap = (float)(qset.size() + mysets->transSets[i].size());
				overlap /= 1 + jaccord_inv;
				//cands[size] = Res(i, overlap);
				cands[size] = Res(i, candRes[i]);
				if (size < ub) {
					size++;
					std::push_heap(cands, cands + size);
				}
				else if (cands[0].inp < cands[size].inp) {
					size++;
					std::push_heap(cands, cands + size);
					std::pop_heap(cands, cands + size);
					size--;
				}
			}

		}
		if (size < ub) ub = size;
		q.res.resize(ub);
		q.time_sift = timer1.elapsed();
		timer1.restart();

		for (int i = 0; i < ub; ++i) q.res[i] = Res(cands[i].id, data->dotProduct2SV(cands[i].id, q.qvec));
		std::partial_sort(q.res.begin(), q.res.begin() + q.k, q.res.end());
		q.res.resize(q.k);

		q.time_verify = timer1.elapsed();
		timer1.restart();

		q.time_total = timer.elapsed();
		return;
	}
};

class myMinHashv3 :public myMinHashBase
{
	const std::string algName = "minHashv3";
	int bucketSize = 0;
	std::vector<std::vector<std::vector<int>>> myIndexes;
	
public:
	myMinHashv3(Preprocess& prep, mips2set& sets, int b_, int r_) :myMinHashBase(prep, sets, b_, r_) {
		bucketSize = len;
		lsh::timer timer;

		std::cout << "\nMINv3 INDEXING..." << std::endl;
		timer.restart();
		buildIndex();
		std::cout << "BUILDING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;
		
	}

	myMinHashv3(myMinHashBase& minBase) :myMinHashBase(minBase) {
		bucketSize = len;
		lsh::timer timer;

		std::cout << "\nMINv3 INDEXING..." << std::endl;
		timer.restart();
		buildIndex();
		std::cout << "BUILDING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

	}

	std::string getAlgName() const {
		return algName;
	}

	int getKeys(int* hs) {
		size_t res = 0;
		for (int j = 0; j < b; ++j) {
			res *= (len % bucketSize);
			//res *= len;
			res %= bucketSize;
			//res = fastrange64(res * len, bucketSize);
			res += hs[j] % bucketSize;
			//res += fastrange64(hs[j], bucketSize);
		}
		return res;
	}

	void buildIndex() {
		//int width = floor(pow(bucketSize, 1.0 / b));
		myIndexes.resize(r);
		for (int j = 0; j < r; ++j) myIndexes[j].resize(bucketSize);
		lsh::progress_display pd(r);

#pragma omp parallel for
		for (int j = 0; j < r; ++j) {
			for (int i = 0; i < n; ++i) {
				for (int k = 0; k < b; ++k) {
					int key = hashvals[i][j * b + k];
					myIndexes[j][key].push_back(i);
				}
			}
			++pd;
		}

	}

	void knnSearch2(Query& q) {
		lsh::timer timer, timer1;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset);
		std::vector<int> qhash(s, len);
		std::vector<bool> unseen(n, true);
		for (auto& id : qset) {
			for (int j = 0; j < s; ++j) if (qhash[j] > hashFuns[j][id]) qhash[j] = hashFuns[j][id];
		}

		std::vector<int> candRes(n, 0);
		size_t cnt = 0;
		for (int j = 0; j < r; ++j) {
			//std::unordered_set<int> cands;
			for (int k = 0; k < b; ++k) {
				int key = qhash[j * b + k];
				for (auto& id : myIndexes[j][key]) candRes[id]++, cnt++;
			}

		}
		q.cnt = cnt;
		q.time_hash = timer1.elapsed();
		timer1.restart();

		//int ub = 2e3;
		int ub = q.k + q.ub;
		if (ub > n) ub = n;
		Res* cands = new Res[ub + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			if (candRes[i] > 0) {
				float jaccord_inv = (float)r / candRes[i];
				//float jaccord_inv = 1.0 / (1 - std::pow(1 - candRes[i] / (float)r, 1.0 / b));
				float overlap = (float)(qset.size() + mysets->transSets[i].size());
				overlap /= 1 + jaccord_inv;
				//cands[size] = Res(i, overlap);
				cands[size] = Res(i, candRes[i]);
				if (size < ub) {
					size++;
					std::push_heap(cands, cands + size);
				}
				else if (cands[0].inp < cands[size].inp) {
					size++;
					std::push_heap(cands, cands + size);
					std::pop_heap(cands, cands + size);
					size--;
				}
			}
		}
		if (size < ub) ub = size;
		q.res.resize(ub);
		q.time_sift = timer1.elapsed();
		timer1.restart();
		q.cost += ub;
		for (int i = 0; i < ub; ++i) q.res[i] = Res(cands[i].id, data->dotProduct2SV(cands[i].id, q.qvec));
		if (q.res.size() > q.k) {
			std::partial_sort(q.res.begin(), q.res.begin() + q.k, q.res.end());
			q.res.resize(q.k);
		}

		q.time_verify = timer1.elapsed();
		timer1.restart();

		q.time_total = timer.elapsed();
		return;   

	}
};

