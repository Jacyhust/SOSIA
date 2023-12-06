#pragma once
#include "def.h"
#include "Preprocess.h"
#include "basis.h"
#include <cmath>
#include <assert.h>
#include <vector>
#include <queue>
#include <cfloat>

struct ivfPair
{
	int id = -1;
	float val = 0.0f;
	ivfPair() = default;
	ivfPair(int id_, float val_) :id(id_), val(val_) {}
};

class Linscan
{
public:
	std::vector<std::vector<ivfPair>> ivfIndex;
	
public:
	size_t n = 0;
	int dim = 0;
	Linscan() = default;

	Linscan(SparseData* data) {
		n = data->n;
		dim = data->dim;
		std::cout << "LINSCAN INDEXING ...\n";
		lsh::timer timer;
		buildIndex(data);
		std::cout << "INDEXING TIME: " << timer.elapsed() << " s.\n";
	}

	void buildIndex(SparseData* data) {
		ivfIndex.resize(dim);
		lsh::progress_display pd(n);
//#pragma omp parallel for
		for (int i = 0; i < n; ++i) {
			SparseVector* pt = data->generateSV(i);
			for (int j = 0; j < pt->nnz; ++j) {
				ivfIndex[pt->indices[j]].emplace_back(i, pt->val[j]);
			}
			++pd;
		}
	}

	Linscan(SparseData* data, int* cands, int len) {
		n = len;
		dim = data->dim;
		buildIndex(data, cands, len);
	}

	void buildIndex(SparseData* data, int* cands, int len) {
		ivfIndex.resize(dim);
		for (int i = 0; i < len; ++i) {
			SparseVector* pt = data->generateSV(cands[i]);
			for (int j = 0; j < pt->nnz; ++j) {
				ivfIndex[pt->indices[j]].emplace_back(cands[i], pt->val[j]);
			}
		}
	}

	Linscan(SparseData* data, Res* cands, int len) {
		n = len;
		dim = data->dim;
		buildIndex(data, cands, len);
	}

	void buildIndex(SparseData* data, Res* cands, int len) {
		ivfIndex.resize(dim);
		for (int i = 0; i < len; ++i) {
			SparseVector* pt = data->generateSV(cands[i].id);
			for (int j = 0; j < pt->nnz; ++j) {
				ivfIndex[pt->indices[j]].emplace_back(i, pt->val[j]);
			}
		}
	}

	void knnSearch(Query& q) {
		lsh::timer timer;

		std::vector<float> score(n, 0.0f);
		for (int j = 0; j < q.queryPoint.nnz; ++j) {
			int ind = q.queryPoint.indices[j];
			float val = q.queryPoint.val[j];
			for (auto& item : ivfIndex[ind]) {
				score[item.id] += val * item.val;
				q.cnt_product++;
			}
		}

		Res* res_PQ= new Res[q.k + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			res_PQ[size] = Res(i, score[i]);
			if (size < q.k) {
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
		int len = size;
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr] = res_PQ[0];
			std::pop_heap(res_PQ, res_PQ + size);
			size--;
			rr--;
		}

		//std::vector<Res> res(q.k + 1);
		q.time_total = timer.elapsed();
	}

	void bruteforce(Query& q, SparseData* data) {
		lsh::timer timer;

		std::vector<float> score(n, 0.0f);
		for (int i = 0; i < n; ++i) {
			SparseVector* sv = data->generateSV(i);
			score[i] = q.queryPoint.dotProduct2SV(sv);
		}

		Res* res_PQ = new Res[q.k + 1];
		int size = 0;
		for (int i = 0; i < n; ++i) {
			res_PQ[size] = Res(i, score[i]);
			if (size < q.k) {
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
		int len = size;
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr] = res_PQ[0];
			std::pop_heap(res_PQ, res_PQ + size);
			size--;
			rr--;
		}

		//std::vector<Res> res(q.k + 1);
		q.time_total = timer.elapsed();
	}
};

