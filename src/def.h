#pragma once

#include <iostream>
#include <vector>
#include <cfloat>
//#include "rd_stdcout.h"
//
//extern rdCoutStream rd;

using DATATYPE = float;
//struct Data
//{
//	// Dimension of data
//	int dim;
//	// Number of data
//	int n;
//	int nq = 0;
//	// Data matrix
//	float** val;
//	float** queries;
//};

//struct SpareseVector;
struct SparseData;

//Define my sparse vector
struct SparseVector
{
	//int dim = 0;
	int nnz = 0;
	int* indices = nullptr;
	float* val = nullptr;

	//friend struct SparseData;

	SparseVector() = default;
	SparseVector(int nnz_,int* ind_,float* val_):
		nnz(nnz_),
		indices(ind_),
		val(val_) {}

	//SparseVector(SparseData* data, int id) :
	//	nnz(data->indptr[id + 1] - data->indptr[id]),
	//	indices(data->indices + data->indptr[id]),
	//	val(data->val + data->indptr[id]) {}

	float getSquareNorm() {
		float res = 0.0f;
		for (int j = 0; j < nnz; ++j) {
			res += val[j] * val[j];
			//data.val[i][indices[j]] = val[j];
		}
		return res;
	}

	float dotProduct2DV(float* dv) {
		float res = 0.0f;
		for (int it = 0; it < nnz; ++it) {
			res += val[it] * dv[indices[it]];
		}
		return res;
	}

	float getMax() {
		float res = -FLT_MAX;
		for (int i = 0; i < nnz; ++i)
			if (res < val[i]) res = val[i];

		return res;
	}

	//compute the inner product to the id-th sparse vector in sd
	//float dotProduct2SV(SparseData* sd, int id) {
	//	float res = 0.0f;

	//	int it1 = 0, it2 = 0;
	//	int* pt1 = sd->indices + sd->indptr[id], * pt2 = indices;
	//	int end1 = sd->indptr[id + 1] - sd->indptr[id], end2 = nnz;
	//	float* val1 = sd->val + sd->indptr[id], * val2 = val;
	//	while (it1 < end1 && it2 < end2) {
	//		if (pt1[it1] == pt2[it2]) res += val1[it1] * val2[it2];
	//		else if (pt1[it1] < pt2[it2]) it1++;
	//		else it2++;
	//	}
	//	return res;
	//}

	//compute the inner product to the sparse vector sv
	float dotProduct2SV(SparseVector* sv) {
		float res = 0.0f;

		int it1 = 0, it2 = 0;
		int* pt1 = sv->indices, * pt2 = indices;
		int end1 = sv->nnz, end2 = nnz;
		float* val1 = sv->val, * val2 = val;
		while (it1 < end1 && it2 < end2) {
			if (pt1[it1] == pt2[it2]) res += val1[it1++] * val2[it2++];
			else if (pt1[it1] < pt2[it2]) it1++;
			else it2++;
		}
		return res;
	}
};

//Define my sparse matrix
struct SparseData
{
	// Dimension of data
	int dim = 0; //ncol
	// Number of data
	size_t n = 0; // nrow
	//int nq = 0;
	// Data matrix
	//float** val;
	size_t nnz = 0;
	size_t* indptr = nullptr; //nrow+1
	int* indices = nullptr; // nnz
	float* val = nullptr; //nnz
	//float** queries;

	SparseData() = default;
	//SparseData(std::string fname){} Read from file, to be implemented in Preprocess.cpp

	SparseVector* generateSV(int id) {
		return new SparseVector(indptr[id + 1] - indptr[id], indices + indptr[id], val + indptr[id]);
	}

	SparseVector generateSVnpt(int id) {
		return SparseVector(indptr[id + 1] - indptr[id], indices + indptr[id], val + indptr[id]);
	}

	float getSquareNorm(int id) {
		int& i = id;
		float res = 0.0f;
		for (int j = indptr[i]; j < indptr[i + 1]; ++j) {
			res += val[j] * val[j];
			//data.val[i][indices[j]] = val[j];
		}
		return res;
	}

	////compute the inner product between the `id`-th point and a sparse vector `sv`
	//float dotProduct2SV(int id, SparseVector& sv) {
	//	float res = 0.0f;
	//	//int it1 = indices[indptr[id]], it2 = sv.indices[0];
	//	//int end1 = indptr[id+1] - indptr[id], end2 = sv.nnz;
	//	//int it1 = indptr[id], it2 = 0;
	//	//int end1 = indptr[id + 1] - indptr[id], end2 = sv.nnz;

	//	int it1 = 0, it2 = 0;
	//	int* pt1 = indices + indptr[id], * pt2 = sv.indices;
	//	int end1 = indptr[id + 1] - indptr[id], end2 = sv.nnz;
	//	float* val1 = val + indptr[id], * val2 = sv.val;
	//	//int it1 = indices[indptr[id]], it2 = sv.indices[0];
	//	//int end1= indices[indptr[id+1]],end
	//	while (it1 < end1 && it2 < end2) {
	//		if (pt1[it1] == pt2[it2]) res += val1[it1] * val2[it2];
	//		else if (pt1[it1] < pt2[it2]) it1++;
	//		else it2++;
	//	}
	//	return res;
	//}

	float getMin() {
		float res = FLT_MAX;
		for (int i = 0; i < nnz; ++i) {
			if (res > val[i]) res = val[i];
		}
		return res;
	}

	float getMax() {
		float res = -FLT_MAX;
		for (int i = 0; i < nnz; ++i) {
			if (res < val[i]) res = val[i];
		}
		return res;
	}

	//compute the inner product between the `id`-th point and a dense vector `dv`
	float dotProduct2DV(int id, float* dv) {
		float res = 0.0f;
		int* pt1 = indices + indptr[id];
		int end1 = indptr[id + 1] - indptr[id];
		float* val1 = val + indptr[id];
		for (int it = 0; it < end1; ++it) {
			res += val1[it] * dv[pt1[it]];
		}
		return res;
	}

	float dotProduct2SV(int id, float* sv) {
		float res = 0.0f;
		int* pt1 = indices + indptr[id];
		int end1 = indptr[id + 1] - indptr[id];
		float* val1 = val + indptr[id];
		for (int it = 0; it < end1; ++it) {
			if (sv[pt1[it]])
				res += val1[it] * sv[pt1[it]];
		}
		return res;
	}

	void showInfo() {
		std::cout << "N=     " << n << "\n";
		std::cout << "dim=   " << dim << "\n\n";
		std::cout << "nnz=   " << nnz << "\n\n";
		std::cout << "Spa.=  " << (float)nnz / (n) << "\n\n";
		std::cout << "Range= " << "[" << getMin() << "," << getMax() << "]" << "\n\n";
	}

	~SparseData() {
		delete[] indices;
		delete[] indptr;
		delete[] val;
	}
};

struct Ben
{
	int N;
	int num;
	int** indice;
	float** innerproduct;
};

struct HashParam
{
	// the value of a in S hash functions
	float** rndAs1;
	// the value of a in S hash functions
	float** rndAs2;
};

