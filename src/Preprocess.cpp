#include "Preprocess.h"
#include "basis.h"
#include <fstream>
#include <assert.h>
#include <random>
#include <iostream>
#include <fstream>
#include <map>
#include <ctime>
#include <sstream>
#include <numeric>
#include<algorithm>
#include <iomanip>

#define CANDIDATES 100
#define E 2.718281746
#define PI 3.1415926
#define MAXSIZE 40960
#define MAX_TOP_k 100


#define min(a,b)            (((a) < (b)) ? (a) : (b))

Preprocess::Preprocess(const std::string& path, const std::string& qpath, const std::string& ben_file_)
{
	lsh::timer timer;
	std::cout << "LOADING DATA..." << std::endl;
	timer.restart();
	data = new SparseData();
	queries = new SparseData();
	load_data(path, data);
	load_data(qpath, queries);
	std::cout << "LOADING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;
	calSquareLen();

	dataFile = path;
	benFile = ben_file_;
	ben_create();
}

//Load Sparse Data
void Preprocess::load_data(const std::string& path, SparseData* sd)
{
	std::ifstream in(path.c_str(), std::ios::binary);
	if (!in) {
		std::cout << "Fail to open the file:\n" << path << "!\n";
		exit(-10086);
	}


	size_t header[3] = {};
	in.read((char*)header, sizeof(header));

	size_t nrow = header[0];
	size_t ncol = header[1];
	size_t nnz = header[2];

	sd->dim = ncol;
	sd->n = nrow;
	sd->nnz = nnz;

	sd->indptr = new size_t[nrow + 1];
	in.read((char*)sd->indptr, sizeof(size_t) * (nrow + 1));

	sd->indices = new int[nnz];
	in.read((char*)sd->indices, sizeof(int) * nnz);

	sd->val = new float[nnz];
	in.read((char*)sd->val, sizeof(float) * nnz);

	std::cout << "Load the sparse data from the file: " << path << "\n";
	sd->showInfo();
	in.close();
}

void Preprocess::calSquareLen()
{
	squareLen = new float[data->n];
	for (int i = 0; i < data->n; ++i) squareLen[i] = data->getSquareNorm(i);
	maxLen = *std::max_element(squareLen, squareLen + data->n);
}

// typedef struct Tuple
// {
// 	int id;
// 	float inp;
// 	bool operator < (const Tuple& rhs) {
// 		return inp < rhs.inp;
// 	}
// }Tuple;
// bool comp(const Tuple& a, const Tuple& b)
// {
// 	return a.inp > b.inp||()
// }

void Preprocess::ben_make()
{
	benchmark.N = queries->n, benchmark.num = MAX_TOP_k;
	benchmark.N = 200;

	if (benchmark.N > queries->n) benchmark.N = queries->n;
	if (benchmark.num > data->n) benchmark.num = data->n;

	benchmark.indice = new int* [benchmark.N];
	benchmark.innerproduct = new float* [benchmark.N];
#pragma omp parallel for
	for (int j = 0; j < benchmark.N; j++){
		benchmark.indice[j] = new int[benchmark.num];
		benchmark.innerproduct[j] = new float[benchmark.num];
	}

	

	lsh::progress_display pd(benchmark.N);
#pragma omp parallel for
	for (int j = 0; j < benchmark.N; j++){
		Res resPair;
		std::vector<Res> dists;
		SparseVector sv=queries->generateSVnpt(j);
		dists.clear();   
		for (int i = 0; i < data->n; i++){
			resPair.id = i;
			resPair.inp = sv.dotProduct2SV(data->generateSV(i));
			dists.emplace_back(resPair);
		}

		sort(dists.begin(), dists.end());
		for (int i = 0; i < benchmark.num; i++){
			benchmark.indice[j][i] = (int)dists[i].id;
			benchmark.innerproduct[j][i] = dists[i].inp;
		}
		++pd;
	}

}



void Preprocess::ben_save()
{
	std::ofstream out(benFile.c_str(), std::ios::binary);
	out.write((char*)&benchmark.N, sizeof(int));
	out.write((char*)&benchmark.num, sizeof(int));

	for (int j = 0; j < benchmark.N; j++) {
		out.write((char*)&benchmark.indice[j][0], sizeof(int) * benchmark.num);
	}

	for (int j = 0; j < benchmark.N; j++) {
		out.write((char*)&benchmark.innerproduct[j][0], sizeof(float) * benchmark.num);
	}

	out.close();
}

void Preprocess::ben_load()
{
	std::ifstream in(benFile.c_str(), std::ios::binary);
	in.read((char*)&benchmark.N, sizeof(int));
	in.read((char*)&benchmark.num, sizeof(int));

	benchmark.indice = new int* [benchmark.N];
	benchmark.innerproduct = new float* [benchmark.N];
	for (int j = 0; j < benchmark.N; j++) {
		benchmark.indice[j] = new int[benchmark.num];
		in.read((char*)&benchmark.indice[j][0], sizeof(int) * benchmark.num);
	}

	for (int j = 0; j < benchmark.N; j++) {
		benchmark.innerproduct[j] = new float[benchmark.num];
		in.read((char*)&benchmark.innerproduct[j][0], sizeof(float) * benchmark.num);
	}
	in.close();
}

void Preprocess::ben_create()
{
	int a_test = data->n + 1;
	lsh::timer timer;
	std::ifstream in(benFile.c_str(), std::ios::binary);
	in.read((char*)&a_test, sizeof(int));
	in.close();
	if (a_test > 0 && a_test <= data->n)
	{
		std::cout << "LOADING BENMARK..." << std::endl;
		timer.restart();
		ben_load();
		std::cout << "LOADING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;
	}
	else
	{
		std::cout << "MAKING BENMARK..." << std::endl;
		timer.restart();
		ben_make();
		std::cout << "MAKING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;

		std::cout << "SAVING BENMARK..." << std::endl;
		timer.restart();
		ben_save();
		std::cout << "SAVING TIME: " << timer.elapsed() << "s." << std::endl << std::endl;
	}


	int qn=100;
	int knn=50;
	// std::cout << "The benchmarks are:\n";
	// for (int i = 0; i < qn; ++i) {
	// 	for (int j = 0; j < knn; ++j)
	// 		std::cout << std::setw(12) << benchmark.indice[i][j] << "\t";
	// 	std::cout << "\n";
	// }

	// for (int i = 0; i < qn; ++i) {
	// 	for (int j = 0; j < knn; ++j)
	// 		std::cout << std::setw(12) << benchmark.innerproduct[i][j] << "\t";
	// 	std::cout << "\n";
	// }
}

Preprocess::~Preprocess()
{
	//clear_2d_array(data.val, data.n);
	//clear_2d_array(data.queries, data.nq);
	clear_2d_array(benchmark.indice, benchmark.N);
	clear_2d_array(benchmark.innerproduct, benchmark.N);
	delete[] squareLen;
}

Partition::Partition(float c_, float c0_, Preprocess& prep)
{
	ratio = (pow(c0_, 4.0f) - 1) / (pow(c0_, 4.0f) - c_);
	makeChunks(prep);
}

Partition::Partition(float c_, Preprocess& prep)
{
	ratio = 0.95;
	float c0_ = 1.5f;
	
	makeChunks(prep);
}



void Partition::makeChunks(Preprocess& prep)
{
	distpairs.clear();
	std::vector<int> bucket;
	Dist_id pair;
	int N_ = prep.data->n;
	int n;
	for (int j = 0; j < N_; j++){
		pair.id = j;
		pair.dist = prep.squareLen[j];
		distpairs.push_back(pair);
	}
	std::sort(distpairs.begin(), distpairs.end());

	num_chunk = 0;
	chunks.resize(N_);
	int j = 0;
	while (j < N_){
		float M = distpairs[j].dist / ratio;
		n = 0;
		bucket.clear();
		while (j < N_){
			if ((distpairs[j].dist > M || n >= MAXSIZE)) {
				break;
			}

			chunks[distpairs[j].id] = num_chunk;
			bucket.push_back(distpairs[j].id);
			j++;
			n++;
		}
		nums.push_back(n);
		MaxLen.push_back(distpairs[(size_t)j - 1].dist);
		EachParti.push_back(bucket);
		bucket.clear();
		num_chunk++;
	}

	newIds.resize(N_, 0);
	std::vector<int> cnts(num_chunk, 0);
	for (int i = 0; i < N_; ++i) {
		newIds[i] = cnts[chunks[i]]++;
	}

	display();
}

void Partition::display()
{
	std::vector<int> n_(num_chunk, 0);
	int N_ = std::accumulate(nums.begin(), nums.end(), 0);
	for (int j = 0; j < N_; j++){
		n_[chunks[j]]++;
	}
	bool f1 = false, f2 = false;
	for (int j = 0; j < num_chunk; j++)
	{
		if (n_[j] != nums[j]) {
			f1 = true;
			break;
		}
	}

	std::cout << "This is the result of partition:"
		<< "\n Blocks       =" << num_chunk
		<< "\n ratio_       =" << sqrt(ratio)
		<< "\n n_pts_       =" << N_ << std::endl;
}

Partition::~Partition(){}


Parameter::Parameter(Preprocess& prep, int L_, int K_, int M_)
{
	N = prep.data->n;
	dim = prep.data->dim;
	L = L_;
	K = K_;
	S = L * K;
	MaxSize = 3;
	KeyLen = 1;
}

Parameter::Parameter(Preprocess& prep, int L_, int K_, int M_, float U_)
{
	N = prep.data->n;
	dim = prep.data->dim;
	M = M_;
	L = L_;
	K = K_;
	S = L * K;
	U = U_;
}

Parameter::Parameter(Preprocess& prep, int L_, int K_, int M_, float U_,float W_)
{
	N = prep.data->n;
	dim = prep.data->dim;
	L = L_;
	K = K_;
	S = L * K;
	U = U_;
	W = W_;
}

Parameter::Parameter(Preprocess& prep, float c_, float S0)
{

	N = prep.data->n;
	dim = prep.data->dim;
	assert(c_ * S0 < 1);
	double pi = atan(1) * 4;
	double p1 = 1 - acos(S0) / pi;
	double p2 = 1 - acos(c_ * S0) / pi;
	K = (int)floor(log(N) / log(1 / p2)) + 1;
	double rho = log(p1) / log(p2);
	L = (int)floor(pow((double)N, rho)) + 1;
	S = L * K;

	
}

inline float normal_pdf0(			// pdf of Guassian(mean, std)
	float x,							// variable
	float u,							// mean
	float sigma)						// standard error
{
	float ret = exp(-(x - u) * (x - u) / (2.0f * sigma * sigma));
	ret /= sigma * sqrt(2.0f * PI);
	return ret;
}

float new_cdf0(						// cdf of N(0, 1) in range [-x, x]
	float x,							// integral border
	float step)							// step increment
{
	float result = 0.0f;
	for (float i = -x; i <= x; i += step) {
		result += step * normal_pdf0(i, 0.0f, 1.0f);
	}
	return result;
}

inline float calc_p0(			// calc probability
	float x)							// x = w / (2.0 * r)
{
	return new_cdf0(x, 0.001f);		// cdf of [-x, x]
}

Parameter::Parameter(Preprocess& prep, float c_)
{
	K = 0;
	L = 0;
	N = prep.data->n;
	dim = prep.data->dim;

	KeyLen = -1;
	MaxSize = -1;
	float w_ = sqrt((8.0f * c_ * c_ * log(c_)) / (c_ * c_ - 1.0f));

	float beta_;
	float delta_;
	float p1_;
	float p2_;
	float para1;
	float para2;
	float para3;
	float eta;
	float alpha_;
	float m, l;

	int n_pts_ = MAXSIZE;
	beta_ = (float)CANDIDATES / n_pts_;
	delta_ = 1.0f / E;

	p1_ = calc_p0(w_ / 2.0f);
	p2_ = calc_p0(w_ / (2.0f * c_));

	para1 = sqrt(log(2.0f / beta_));
	para2 = sqrt(log(1.0f / delta_));
	para3 = 2.0f * (p1_ - p2_) * (p1_ - p2_);
	eta = para1 / para2;

	alpha_ = (eta * p1_ + p2_) / (1.0f + eta);
	m = (para1 + para2) * (para1 + para2) / para3;
	this->S = (int)ceil(m);
	M = 1;
	W = 1;
}


bool Parameter::operator = (const Parameter& rhs)
{
	bool flag = 1;
	flag *= (L == rhs.L);
	flag *= (K = rhs.K);
	return flag;
}

Parameter::~Parameter()
{}