#pragma once
#include "method.h"

class Sosia :public myMinHashBase
{
	const std::string algName = "sosia";
	int bucketSize = 0;
	std::vector<std::vector<std::vector<int>>> myIndexes;

	std::vector<Res> sortedNorms;

public:
	Sosia(Preprocess& prep, mips2set& sets, int b_, int r_) :myMinHashBase(prep, sets, b_, r_) {
		bucketSize = len;
		lsh::timer timer;

		std::cout << "\nSOSIA INDEXING..." << std::endl;
		timer.restart();
		sortNorm(sets);
		buildIndex();
		std::cout << "BUILDING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

	}

	Sosia(myMinHashBase& minBase) :myMinHashBase(minBase) {
		bucketSize = len;
		lsh::timer timer;

		std::cout << "\nMINv3 INDEXING..." << std::endl;
		timer.restart();
		sortNorm(*mysets);
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

	//void sortNorm(Preprocess& prep) {
	//	sortedNorms.resize(n);
	//	auto data = prep.data;
	//	for (int i = 0; i < n; ++i) {
	//		sortedNorms[i] = Res(i, data->generateSV(i)->getSquareNorm());
	//	}
	//}

	void sortNorm(mips2set& sets) {
		sortedNorms.resize(n);
		//auto data = set;
		for (int i = 0; i < n; ++i) {
			sortedNorms[i] = Res(i, sets.transSets[i].size());
		}
		sort(sortedNorms.begin(), sortedNorms.end());
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
					int id = sortedNorms[i].id;
					int key = hashvals[id][j * b + k];
					myIndexes[j][key].push_back(i);
					//myIndexes[j][key].push_back(sortedNorms[i].id);
				}
			}

			//for (int i = 0; i < bucketSize; ++i) {
			//	sort(myIndexes[j][i].begin(), myIndexes[j][i].end());
			//}

			++pd;
		}
	}

	void knnSearch2(Query& q) {
		lsh::timer timer, timer1;
		std::vector<int> qset;
		mips2set::getQuerySet(q, qset); // the reason why it is random
		std::vector<int> qhash(s, len);
		std::vector<bool> unseen(n, true);
		for (auto& id : qset) {
			for (int j = 0; j < s; ++j) if (qhash[j] > hashFuns[j][id]) qhash[j] = hashFuns[j][id];
		}

		std::vector<int> candRes(n, 0);
		size_t cnt = 0;
		int qm=q.qm;
		if(q.qm<0) qm=r;
		for (int j = 0; j < qm; ++j) {
			//std::unordered_set<int> cands;
				int key = qhash[j * b];
				for (auto& id : myIndexes[j][key]) candRes[id]++, cnt++;
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
				//if((qset.size() + sortedNorms[i].inp)*candRes[i]>ptau*m*)
				float jaccord_inv = (float)qm / candRes[i];
				//float jaccord_inv = 1.0 / (1 - std::pow(1 - candRes[i] / (float)r, 1.0 / b));
				float overlap = qset.size() + sortedNorms[i].inp;
				float L = overlap;
				overlap /= 1 + jaccord_inv;
				cands[size] = Res(sortedNorms[i].id, overlap);
				//cands[size] = Res(sortedNorms[i].id, candRes[i] * L);
				//cands[size] = Res(sortedNorms[i].id, candRes[i]);
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