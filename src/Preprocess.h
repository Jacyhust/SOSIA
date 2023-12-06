#pragma once
#include "def.h"
#include <cmath>
#include <assert.h>
#include <string>
#include <vector>
class Preprocess
{
public:
	SparseData* data = nullptr;
	SparseData* queries = nullptr;
	float* squareLen = nullptr;
	Ben benchmark;
	float maxLen;
	std::string dataFile;
	std::string benFile;
public:
	Preprocess(const std::string& path, const std::string& qpath, const std::string& ben_file_);
	void load_data(const std::string& path, SparseData* sd);
	//void load_query(const std::string& path);
	void calSquareLen();
	void ben_make();
	void ben_save();
	void ben_load();
	void ben_create();
	~Preprocess();
};

struct Dist_id
{
	int id;
	float dist;
	bool operator < (const Dist_id& rhs) {
		return dist < rhs.dist;
	}
};

class Partition
{
private:
	float ratio;
	void makeChunks(Preprocess& prep);
public:
	int num_chunk;
	std::vector<float> MaxLen;
	//The chunk where each point belongs
	std::vector<int> chunks;
	//The data size of each chunks
	std::vector<int> nums;
	//the new id for n points
	std::vector<int> newIds;    
	//The buckets by parti;
	std::vector<std::vector<int>> EachParti;
	std::vector<Dist_id> distpairs;
	void display();

	Partition(float c_, Preprocess& prep);

	Partition(float c_, float c0_, Preprocess& prep);
	//Partition() {}
	~Partition();
};

class Parameter //N,dim,S, L, K, M, W;
{
public:
	int N;
	int dim;
	// Number of hash functions
	int S;
	//#L Tables; 
	int L;
	// Dimension of the hash table
	int K;
	//
	int MaxSize;
	//
	int KeyLen;

	int M = 1;

	int W = 0;

	float U;
	Parameter(Preprocess& prep, int L_, int K_, int M);
	Parameter(Preprocess& prep, int L_, int K_, int M_, float U_);
	Parameter(Preprocess& prep, int L_, int K_, int M_, float U_, float W_);
	Parameter(Preprocess& prep, float c_, float S0);
	Parameter(Preprocess& prep, float c0_);
	bool operator = (const Parameter& rhs);
	~Parameter();
};


struct Res//the result of knns
{
	float inp;
	int id;
	Res() = default;
	Res(int id_, float inp_) :id(id_), inp(inp_) {}
	bool operator < (const Res& rhs) const {
		return inp > rhs.inp||(inp == rhs.inp && id < rhs.id);
	}
};

struct ResInt//the result of knns
{
	int inp = 0;
	int id = 0;
	ResInt() = default;
	ResInt(int id_, float inp_) :id(id_), inp(inp_) {}
	bool operator < (const ResInt& rhs) const {
		return inp > rhs.inp;
	}

	bool operator > (const ResInt& rhs) const {
		return inp < rhs.inp;
	}
};

struct Query
{
	// the parameter "c" in "c-ANN"
	float c;
	SparseVector queryPoint;
	float* qvec = nullptr;
	// k-NN
	int k;
	// Indice of query point in query set.
	int qid;
	//
	int ub;
	//
	size_t cnt = 0;
	size_t cnt_product = 0;
	//
	int cost = 0;
	//
	float time_total = 0;
	//
	float time_hash = 0;
	//
	float time_sift = 0;

	float time_verify = 0;
	// query result: <indice of ANN, distance of ANN>
	std::vector<Res> res;

public:
	int qm=-1;
	Query() = default;
	Query(int id, float c_, int k_, Preprocess& prep, int ub_) {
		qid = id;
		c = c_;
		k = k_;
		ub = ub_;
		queryPoint = prep.queries->generateSVnpt(qid);
		qvec = new float[prep.data->dim];
		for (int i = 0; i < prep.data->dim; ++i) {
			qvec[i] = 0;
		}
		for (int i = 0; i < queryPoint.nnz; ++i) {
			qvec[queryPoint.indices[i]] = queryPoint.val[i];
		}

	}

	~Query() {
		delete[] qvec;
	}
};