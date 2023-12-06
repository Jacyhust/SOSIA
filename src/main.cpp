

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <string.h>
#include <cstring>
#include <chrono>
#include <thread>
//#include <pthread.h>
#include "Preprocess.h"
#include "mf_alsh.h"
#include "basis.h"
#include "alg.h"
#include "expe.h"
#include "rd_stdcout.h"


extern std::string data_fold, index_fold;
extern std::string data_fold1, data_fold2;
int mips2set::l = 10;
std::mt19937 mips2set::rng(int(std::time(0)));
std::uniform_real_distribution<float> mips2set::ur(0, 1);

std::string getNameByDate();

//Load VBMW Data
using dt = uint32_t;
void loadVBMW_dataset(const std::string& path)
{
	dt* res = nullptr;
	std::ifstream in(path.c_str(), std::ios::binary);
	if (!in) {
		std::cout << "Fail to open the file:\n" << path << "!\n";
		exit(-10086);
	}
	//dt n=0;
	//in.read((char*)&n, sizeof(dt));
	//std::cout << "n= [" << n << "]: " << std::endl;
	dt* tmp = new dt[100000];
	int cnt = 0;
	std::vector<dt> lens;
	while (!in.eof()) {
		int len = 0;
		in.read((char*)&len, sizeof(dt));
		in.read((char*)tmp, sizeof(dt) * (len));
		lens.push_back(len);
		if (cnt == 0) {
			std::cout << cnt++ << "-th posting list [" << len << "]: " << std::endl;
			for (int i = 0; i < len; ++i) {
				std::cout << tmp[i] << " ";
			}
			std::cout << std::endl;
		}

		cnt++;
		
	}
	std::cout << "cnt= [" << cnt << "]: " << std::endl;
	return;
}

int main(int argc, char const* argv[])
{
	rdCoutStream rd(getNameByDate()); //Redirect the cout to a file and output it to the console at the same time (Create at Nov 13, 2023)

	std::string dataset = "base_small";
	std::string qname= "base_queries";

	int datasetID = -1;
	
	//dataset = "test";
	//qname = "test_queries";
	
	int L = 5;
	int K = 12;
	int k = 50;
	float c = 0.5;
	int l = 20;
	int ub = 10000;
	int m = 150;
	int b = 1;
	float F = 2.0f;
	int sin_m = 30, sin_h = 2;
	int mode = 0;
	const int num_th = thread::hardware_concurrency();
	
	if (argc > 1) datasetID = std::stoi(argv[1]);

	switch(datasetID){
	case 0:
		dataset = "base_1M";
		//l = 40;
		ub = 5000;
		m=150;
		l=40;
		qname= "base_queries";
		break;

	case 1:
		dataset = "base_full";
		l = 40;
		ub = 10000;
		m=150;
		qname= "base_queries";
		break;

	case 2:
		dataset = "rand_1M";
		qname= "rand";
		ub = 5000;
		m=150;
		l=40;
		break;

	case 3:
		dataset = "rand_10M";
		qname= "rand";
		l = 40;
		ub = 10000;
		m=150;
		break;

	default:
		break;
	}


	// switch(datasetID){
	// case 0:
	// 	dataset = "base_full";
	// 	l = 40;
	// 	ub = 10000;
	// 	m=150;
	// 	qname= "base_queries";
	// 	break;

	// case 1:
	// 	dataset = "base_1M";
	// 	//l = 40;
	// 	ub = 5000;
	// 	qname= "base_queries";
	// 	break;

	// case 2:
	// 	dataset = "tfidf";
	// 	qname= "tfidf_queries";
	// 	break;
	// case 3:
	// 	dataset = "tfidf1M";
	// 	qname= "tfidf_queries1M";
	// 	break;
	// case 4:
	// 	dataset = "tfidfm1M";
	// 	qname= "tfidfm_queries1M";
	// 	break;
	// case 5:
	// 	dataset = "uni1M";
	// 	qname= "uni1M";
	// 	break;
	// case 6:
	// 	dataset = "uni_full";
	// 	qname= "unim_full";
	// 	break;
	// case 7:
	// 	dataset = "rand_1M";
	// 	qname= "rand";
	// 	ub=5000;
	// 	// dataset = "spladeR_full";
	// 	// qname= "spladeR_full";
	// 	break;
	// case 8:
	// 	dataset = "rand_10M";
	// 	qname= "rand";
	// 	l = 40;
	// 	ub = 10000;
	// 	m=150;
	// 	break;
	// default:
	// 	break;
	// }

	if (argc > 2) l = std::atoi(argv[2]);
	if (argc > 3) mode = std::atoi(argv[3]);

	std::string argvStr[5];
	argvStr[1] = "RawDataSet/sparse/" + dataset + ".csr";
	argvStr[2] = (dataset + ".index");
	argvStr[3] = (dataset + "_all.ben");
	argvStr[4] = "RawDataSet/sparse/" + qname + ".dev.csr";

	std::cout << "Using SOSIA for Sparse Dataset " << argvStr[1] << std::endl;
	std::cout << "l=       " << l << std::endl;
	Preprocess prep(data_fold1 + (argvStr[1]), data_fold1 + (argvStr[4]), data_fold2 + (argvStr[3]));

	std::vector<resOutput> res;

	RecallTime(prep, c, k, m, ub, F, l, sin_m, sin_h, res);
	varyk(prep, c, k, m, ub, F, l, sin_m, sin_h, res);
	//runAll(prep, c, k, m, ub, F, l, sin_m, sin_h, res);

	// varyL(prep, c, k, m, ub, F, l, sin_m, sin_h, res);
	// varyM(prep, c, k, m, ub, F, l, sin_m, sin_h, res);





	saveAndShow(c, k, dataset, res);
	return 0;
}
