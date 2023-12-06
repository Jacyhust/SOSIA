#pragma once
#pragma once
#include "linscan.h"
class Sinnamon : public Linscan
{
	//std::vector<float> ubTerms;
	//std::vector<int> wandIters;
	//std::vector<int> docIters;
	//float thred = 0.0f;
	//int pivot = -1;
	SparseData* data = nullptr;
	int m = 0;
	int h = 0;
	std::vector<std::vector<float>> ubs;
	std::vector<std::vector<float>> lbs;
	std::vector<std::vector<int>> hashFuns;
	
	//int** hashvals = nullptr;
public:
	float T=0.0f;
	Sinnamon(SparseData* data, int m_, int h_, float T_) : Linscan(data) {
		this->data = data;
		m = m_;
		h = h_;
		T = T_;
		init();
	}

	Sinnamon(Linscan& lin, SparseData* data_, int m_, int h_, float T_) : Linscan(lin) {
		this->data = data_;
		m = m_;
		h = h_;
		T = T_;
		init();
	}

	void init() {
		lsh::timer timer;
		std::cout << "SETTING HASH PARAMETER..." << std::endl;
		timer.restart();
		setHash();
		std::cout << "SETTING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

		std::cout << "COMPUTING HASH..." << std::endl;
		timer.restart();
		getHash();
		std::cout << "COMPUTING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;
	}

	void setHash() {
		hashFuns.resize(h);
		//lsh::progress_display pd(h);
#pragma omp parallel for
		for (int i = 0; i < h; ++i) {
			hashFuns[i].resize(dim);
			unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
			for (int j = 0; j < dim; ++j) hashFuns[i][j] = j;
			std::shuffle(hashFuns[i].begin(), hashFuns[i].begin() + dim, std::default_random_engine(seed));
			for (int j = 0; j < dim; ++j) hashFuns[i][j] %= m;
			//++pd;
		}
	}

	void getHash() {
		ubs.resize(n, std::vector<float>(m, 0));//0? or FLT_MAX?
		lbs.resize(n, std::vector<float>(m, FLT_MAX));
		//for (int i = 0; i < n; ++i)hashvals[i] = new int[s];

		lsh::progress_display pd(n);
#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			SparseVector x = data->generateSVnpt(i);
			for (int j = 0; j < h; ++j) {
				for (int k = 0; k < x.nnz; ++k) {
					int ind = x.indices[k];
					float val = x.val[k];
  					
					if (val > ubs[i][hashFuns[j][ind]]) ubs[i][hashFuns[j][ind]] = val;
					if (val < ubs[i][hashFuns[j][ind]]) lbs[i][hashFuns[j][ind]] = val;
				}
			}
			++pd;
		}
	}

	struct sparsePair
	{
		float val;
		int ind;
		sparsePair() = default;
		sparsePair(int id_, float inp_) :ind(id_), val(inp_) {}
		bool operator < (const sparsePair& rhs) const {
			return val > rhs.val;
		}
	};

	void knnSearch(Query& q) {
		lsh::timer timer, timer1;
		std::vector<float> scores(n, 0.0f);
		std::vector<sparsePair> sp(q.queryPoint.nnz);
		for (int j = 0; j < q.queryPoint.nnz; ++j) {
			int ind = q.queryPoint.indices[j];
			float val = q.queryPoint.val[j];
			sp[j] = sparsePair(ind, val);
		}
		std::sort(sp.begin(), sp.end());

		for (auto& x : sp) {
			for (auto& item : ivfIndex[x.ind]) {
				float ubt = FLT_MAX;
				for (int i = 0; i < h; ++i) {
					if (ubt > ubs[item.id][hashFuns[i][x.ind]]) ubt = ubs[item.id][hashFuns[i][x.ind]];
				}
				scores[item.id] += ubt * x.val;
				if(timer.elapsed()>T) break;
			}
			if(timer.elapsed()>T) break;
		}
		Res* res_PQ = new Res[q.ub + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			res_PQ[size] = Res(i, scores[i]);
			if (size < q.ub) {
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
		if (q.ub > n) q.ub = n;
		q.res.resize(q.ub);
		for (int i = 0; i < q.ub; ++i) {
			q.res[i] = Res(res_PQ[i].id, data->dotProduct2SV(res_PQ[i].id, q.qvec));
		}
		std::partial_sort(q.res.begin(), q.res.begin() + q.k, q.res.end());
		q.res.resize(q.k);

		//q.time_verify = timer1.elapsed();
		//timer1.restart();

		q.time_total = timer.elapsed();
		return;

	}
};