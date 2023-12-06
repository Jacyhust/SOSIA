#pragma once
#include "linscan.h"

struct docTermPair {
	int docid = -1;
	int termid = -1;

	docTermPair() = default;
	docTermPair(int did, int tid) :docid(did), termid(tid) {}
	bool operator < (const docTermPair& rhs) const {
		return docid > rhs.docid;
	}

	bool operator > (const docTermPair& rhs) const {
		return docid < rhs.docid;
	}
};

struct mywandheap {
private:
	int l = 4096;
	docTermPair container[4096];
	//docTermPair container[l];
	int size = 0;

public:
	void emplace(int did, int tid) {
		if (size < l) {
			container[size++] = docTermPair(did, tid);
			std::push_heap(container, container + size);
		}

		else std::cerr << "Out of Index in wandheap!\n";
	}

	bool empty() { return size == 0; }

	docTermPair top() {
		if (size) return container[0];
		else std::cerr << "Out of Index in wandheap!\n";
	}

	void pop() {
		if (size) std::pop_heap(container, container + size), size--;
		else std::cerr << "Out of Index in wandheap!\n";
	}
};

class wand : public Linscan
{
	
	std::vector<float> ubTerms;
	std::vector<int> wandIters;
	std::vector<int> docIters;
	float thred = 0.0f;
	int pivot = -1;
	SparseData* data = nullptr;
	
	//using wandheap = std::priority_queue<docTermPair, std::vector<docTermPair>>;
	//using wandheap = std::priority_queue<docTermPair, std::vector<docTermPair>, std::greater<docTermPair>>;

	//
	//template <class wandheap>
	//inline void wandIterate(wandheap& wandheap, int ind);

	//template <class wandheap>
	//void wandIterate(wandheap& wandheap, int ind) {
	//	if (!wandheap.empty())wandheap.pop();
	//	if (wandIters[ind] < ivfIndex[ind].size()) {
	//		wandheap.emplace(ivfIndex[ind][wandIters[ind]++].id, ind);
	//	}
	//}

public:
	float F = 1.0f;
	wand(SparseData* data) : Linscan(data) {
		this->data = data;
		ubTerms.resize(dim, -1.0f);
		for (int i = 0; i < dim; ++i) {
			for (auto& term : ivfIndex[i]) {
				if (term.val > ubTerms[i]) ubTerms[i] = term.val;
			}
		}
	}

	wand(Linscan& lin, SparseData* data_, float F_=1.0f) : Linscan(lin) {
		F = F_;
		this->data = data_;
		ubTerms.resize(dim, -1.0f);
		for (int i = 0; i < dim; ++i) {
			for (auto& term : ivfIndex[i]) {
				if (term.val > ubTerms[i]) ubTerms[i] = term.val;
			}
		}
	}

	void knnSearch0(Query& q) {
		lsh::timer timer;
		std::priority_queue<Res> resheap;
		thred = 0.0f;
		pivot = -1;

		// Bug2023.11.16: resize(dim,0) will not set 0 for the required items
		wandIters.resize(dim, 0);
		docIters.resize(dim, 0);
		for (int i = 0; i < dim; ++i) {
			wandIters[i] = 0;
			docIters[i] = 0;
		}


		std::vector<float> ubDocs(n, 0.0f);
		std::vector<float> ubq(ubTerms);
		//sort(terms,posting)
		//using wandheap = mywandheap;
		mywandheap wandheap;
		//wandheap.reserve(100);
		//init
		//if (q.qid == 3)
		//	int a = 1;
		for (int j = 0; j < q.queryPoint.nnz; ++j) {
			int ind = q.queryPoint.indices[j];
			float val = q.queryPoint.val[j];
			ubq[ind] *= val;
			if (wandIters[ind] < ivfIndex[ind].size()) {
				wandheap.emplace(ivfIndex[ind][wandIters[ind]++].id, ind);
			}
		}

		//if (q.qid == 3)
		//	int a = 1;

		while (!wandheap.empty()) {
			auto top = wandheap.top();
			//wandIterate(wandheap, top.termid);

			auto& ind = top.termid;
			if (!wandheap.empty())wandheap.pop();
			if (wandIters[ind] < ivfIndex[ind].size()) {
				wandheap.emplace(ivfIndex[ind][wandIters[ind]++].id, ind);
			}

			if (top.docid <= pivot) continue;

			ubDocs[top.docid] += ubq[top.termid];
			if (ubDocs[top.docid] > thred * F) {
				pivot = top.docid;
				//updateRes(resheap, q);
				updateResFast(resheap, q);
			}

		}


		int len = resheap.size();
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr--] = resheap.top();
			resheap.pop();
		}

		q.time_total = timer.elapsed();
	}

	void knnSearch1(Query& q) {
		lsh::timer timer;
		std::priority_queue<Res> resheap;
		thred = 0.0f;
		pivot = -1;

		// Bug2023.11.16: resize(dim,0) will not set 0 for the required items
		wandIters.resize(dim, 0);
		docIters.resize(dim, 0);
		for (int i = 0; i < dim; ++i) {
			wandIters[i] = 0;
			docIters[i] = 0;
		}

		std::vector<float> ubDocs(n, 0.0f);
		std::vector<float> ubq(ubTerms);
		//sort(terms,posting)
		//using wandheap = std::priority_queue<docTermPair, std::vector<docTermPair>>;
		//wandheap wandheap;
		std::priority_queue<docTermPair, std::vector<docTermPair>> wandheap;
		//init
		//if (q.qid == 3)
		//	int a = 1;
		for (int j = 0; j < q.queryPoint.nnz; ++j) {
			int ind = q.queryPoint.indices[j];
			float val = q.queryPoint.val[j];
			ubq[ind] *= val;
			if (wandIters[ind] < ivfIndex[ind].size()) {
				wandheap.emplace(ivfIndex[ind][wandIters[ind]++].id, ind);
			}
		}

		//if (q.qid == 3)
		//	int a = 1;

		while (!wandheap.empty()) {
			auto top = wandheap.top();
			//wandIterate(wandheap, top.termid);

			auto& ind = top.termid;
			if (!wandheap.empty())wandheap.pop();
			if (wandIters[ind] < ivfIndex[ind].size()) {
				wandheap.emplace(ivfIndex[ind][wandIters[ind]++].id, ind);
			}

			if (top.docid <= pivot) continue;

			ubDocs[top.docid] += ubq[top.termid];
			if (ubDocs[top.docid] > thred * F) {
				pivot = top.docid;
				//updateRes(resheap, q);
				updateResFast(resheap, q);
			}

		}


		int len = resheap.size();
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr--] = resheap.top();
			resheap.pop();
		}

		q.time_total = timer.elapsed();
	}

	

	inline void updateRes(std::priority_queue<Res>& resheap, Query& q) {
		q.cost++;
		float score = 0.0f;
		for (int j = 0; j < q.queryPoint.nnz; ++j) {
			int ind = q.queryPoint.indices[j];
			while (docIters[ind] < ivfIndex[ind].size() && ivfIndex[ind][docIters[ind]].id < pivot) docIters[ind]++;
			if (docIters[ind] < ivfIndex[ind].size() && ivfIndex[ind][docIters[ind]].id == pivot)
				score += ivfIndex[ind][docIters[ind]++].val * q.queryPoint.val[j], q.cnt_product++;
		}

		resheap.emplace(pivot, score);
		if (resheap.size() > q.k) resheap.pop();
		if (resheap.size() == q.k) thred = resheap.top().inp;
	}

	void knnSearch2(Query& q) {
		lsh::timer timer;
		std::priority_queue<Res> resheap;
		thred = 0.0f;
		pivot = -1;

		// Bug2023.11.16: resize(dim,0) will not set 0 for the required items
		wandIters.resize(dim, 0);
		docIters.resize(dim, 0);
		for (int i = 0; i < dim; ++i) {
			wandIters[i] = 0;
			docIters[i] = 0;
		}
		std::vector<float> ubDocs(n, 0.0f);
		std::vector<float> ubq(ubTerms);
		for (int j = 0; j < q.queryPoint.nnz; ++j) {
			int ind = q.queryPoint.indices[j];
			float val = q.queryPoint.val[j];
			ubq[ind] *= val;
			for (auto& item : ivfIndex[ind]) {
				ubDocs[item.id] += ubq[ind];
			}
		}

		for (int i = 0; i < n; ++i) {
			if (ubDocs[i] > thred * F) {
				pivot = i;
				//updateResFast(resheap, q);
				updateRes(resheap, q);
			}
		}

		int len = resheap.size();
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr--] = resheap.top();
			resheap.pop();
		}

		q.time_total = timer.elapsed();
	}

	void knnSearch3(Query& q) {
		lsh::timer timer;
		std::priority_queue<Res> resheap;
		thred = 0.0f;
		pivot = -1;

		// Bug2023.11.16: resize(dim,0) will not set 0 for the required items
		wandIters.resize(dim, 0);
		docIters.resize(dim, 0);
		for (int i = 0; i < dim; ++i) {
			wandIters[i] = 0;
			docIters[i] = 0;
		}
		std::vector<float> ubDocs(n, 0.0f);
		std::vector<float> ubq(ubTerms);
		for (int j = 0; j < q.queryPoint.nnz; ++j) {
			int ind = q.queryPoint.indices[j];
			float val = q.queryPoint.val[j];
			ubq[ind] *= val;
			for (auto& item : ivfIndex[ind]) {
				ubDocs[item.id] += ubq[ind];
			}
		}

		for (int i = 0; i < n; ++i) {
			if (ubDocs[i] > thred * F) {
				pivot = i;
				updateResFast(resheap, q);
				//updateRes(resheap, q);
			}
		}

		int len = resheap.size();
		q.res.resize(len);
		int rr = len - 1;
		while (rr >= 0) {
			q.res[rr--] = resheap.top();
			resheap.pop();
		}

		q.time_total = timer.elapsed();
	}


	inline void updateResFast(std::priority_queue<Res>& resheap, Query& q) {
		q.cost++;
		float score = data->dotProduct2SV(pivot, q.qvec);

		resheap.emplace(pivot, score);
		if (resheap.size() > q.k) resheap.pop();
		if (resheap.size() == q.k) thred = resheap.top().inp;
	}
};