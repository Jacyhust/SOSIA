#pragma once
#include "def.h"
#include "Preprocess.h"
#include <cmath>
#include <assert.h>
#include <vector>
#include <queue>
#include <cfloat>

namespace fargo
{
	class Hash
	{
	private:
		std::string index_file;
		
	public:
		int N;
		int dim;
		// Number of hash functions
		int S;
		//#L Tables; 
		int L;
		// Dimension of the hash table
		int K;

		float** hashval;
		Partition parti;
		HashParam hashpar;
		std::vector<int>*** myIndexes;

		float tmin;
		float tstep;
		float smin;
		float sstep;
		int rows;
		int cols;
		float** phi;

		void load_funtable(const std::string& file);
	public:
		Hash(Preprocess& prep_, Parameter& param_, const std::string& file, Partition& part_, const std::string& funtable);
		void SetHash();
		void GetHash(Preprocess& prep);
		void GetTables(Preprocess& prep);
		bool IsBuilt(const std::string& file);
		~Hash();
	};

	struct hash_pair
	{
		float val;
		int bias;

		bool operator < (const hash_pair& rhs) const {
			return val < rhs.val;
		}
	};

	struct indice_pair
	{
		int key;
		int end;
		int table_id;
		float score;
		bool operator < (const indice_pair& rhs) const {
			//return score > rhs.score;
			return score > rhs.score
				//|| (score == rhs.score && key > rhs.key)
				//|| (score == rhs.score && key == rhs.key && table_id > rhs.table_id)
				;
		}
	};

	struct mp_pair
	{
		//int key;
		int end;
		//int table_id;
		float score = FLT_MAX;
	};

	struct gmp_pair
	{
		int key;
		//int end;
		int table_id;
		float score;
		bool operator < (const gmp_pair& rhs) const {
			//return score > rhs.score;
			return score < rhs.score
				//|| (score == rhs.score && key > rhs.key)
				//|| (score == rhs.score && key == rhs.key && table_id > rhs.table_id)
				;
		}

		//bool operator > (const gmp_pair& rhs) const {
		//	//return score > rhs.score;
		//	return score < rhs.score
		//		//|| (score == rhs.score && key > rhs.key)
		//		//|| (score == rhs.score && key == rhs.key && table_id > rhs.table_id)
		//		;
		//}
	};

	class Query
	{
	private:
		// the parameter "c" in "c-ANN"
		float c;
		//which chunk is accessed
		int chunks;
		int UB = 100;
		
		//float* query_point;

		SparseVector queryPoint;

		// the hash value of query point
		float* hashval;

		std::vector<std::vector<hash_pair>> weigh;
		std::vector<float> total_score;

		//float** mydata;
		SparseData* mydata = nullptr;

		int dim;

		float inp_LB;
		// Set of points sifted
		//std::vector<int> candidate;

		std::vector<int> keys;

		void shift(indice_pair& ip0, indice_pair& res);
		void expand(indice_pair& ip0, indice_pair& res);

		std::priority_queue<indice_pair> global_min;
		indice_pair* ProbingSequence;
		int SequenceLen;
		//int tid = -1;
		float varphi(float x, float theta, Hash& hash);
	public:
		// k-NN
		int k;
		// Indice of query point in query set.
		int qid;
		//
		float norm;
		//
		int cost = 0;
		//cost of each partition
		std::vector<int> costs;
		size_t cnt = 0;
		size_t cnt_product = 0;
		//
		float time_total = 0;
		//
		float time_hash = 0;
		//
		float time_sift = 0;

		float time_verify = 0;
		// query result:<indice of ANN,distance of ANN>
		std::vector<Res> res;

		void cal_hash(Hash& hash, Preprocess& prep);
		void sift(Hash& hash, Preprocess& prep);
		void knn(std::priority_queue<Res>& res_PQ,
			Hash& hash, Preprocess& prep,
			std::vector<int>** table,
			std::vector<bool>& flag, int& num_res);
		void siftF(Hash& hash, Preprocess& prep);
		void knnF(Res* res_PQ, Hash& hash, Preprocess& prep, std::vector<int>** table, std::vector<bool>& flag_, int& size);
	public:
		Query(int id, float c_, int k_, Hash& hash, Preprocess& prep, int ub_);

		~Query();
	};
}


