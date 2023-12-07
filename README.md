
# The Source Code for SOSIA (Submitted to ICDE 2024)
-----------------------------------------------------------------------------------------------------------------
## Introduction
This is a source code for the algorithm described in the paper **[Efficient Approximate Maximum Inner Product Search over Sparse Vectors (Submitted to ICDE 2024)]**. We call it as **sos** project.
## Compilation
The **sos** project is written by **C++** (under `C++17` standard) and is simple and easy to use. It can be complied by **g++** in **Linux** or **MSVC** in **Windows**. To completely support `C++17` standard, the **g++** version is suggested to be **version 9** and **MSVC** version is suggested to be at least **MSVC 19.15 (Visual Studio 2017 15.8)**.

### Installation
#### Windows
We can use **Visual Studio 2019** (Other version of Visual Studio should also work but remains untested) to build the project with importing all the files in the directory `./src/`.

#### Linux
```bash
make
```
The excutable file is then in the main directory, called as `sos`

### Command Usage

-------------------------------------------------------------------
> ./sos mode datasetName queriesName
-------------------------------------------------------------------
(the first parameter specifies the procedure be executed and change)

FOR EXAMPLE, YOU CAN RUN THE FOLLOWING CODE IN COMMAND LINE AFTER BUILD ALL THE TOOLS:

```bash
cd .
make
./sos 0 base_small queries
```

#### Parameter explanation

- mode         : 0-4, an integer, the running mode. mode `0` run the `sos` in the default parameters; mode `1` run the `sos` by varying the base `l`; mode `2` run the `sos` by varying the number of minHash functions `m`; mode `3` run the `sos` by varying the `k` of top-k results; mode `4` run the `sos` for obtaining the recall-time curves.
- datasetName  : the dataset file name
- queries      : the query file name
-------------------------------------------------------------------

## Dataset

In our project, the format of the input file (such as `base_small.csr`) is a binary file but not a text file, because binary file has many advantages. The binary file is organized as the following format:

>{nrow(int64)} {ncol(int64)}{nnz(int64)} {indptr\[0:(nrow)\](int64)} {indices\[0:nnz-1\](int32)} {data\[0:nnz-1\](float32)}

### Parameter explanation

- nrow        : a positive integer, the cardinality of the dataset
- ncol        : a positive integer, the dimensionality of the dataset
- nnz         : a positive integer, the number of non-zero values in the datasets
- indptr      : an array with size `nrow+1`. `indptr[i]` stores the number of non-zero values of the first `i-1` points in the dataset. `indptr[0]=0`
- indices     : an array with size `nnz`. `indices` stores all the non-zero indices in the dataset by the order of points
- data        : an array with size `nnz`. `indices` stores all the non-zero values in the dataset by the order of points

FOR EXAMPLE, for the dataset $D=\{x_0,x_1\}$ where $x_0=(0.1,0.2,0)$ and $x_1=(0.3,0,0.5)$, `nrow=2`, `ncol=3` and `nnz=4`. `indptr=[0,2,4]`, `indices=[0,1,0,2]` and `data=[0.1,0.2,0.3,0.5]`. So, the input file should be as following:

> 2(int64)3(int64)4(int64)0(int64)...0.5(float32)

The input file is read via `C++` as following. You can understand the data format via this `C++` function.
```c++
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
```

A sample dataset `base_small` and `queries.dev.csr` have been put in the directory `./datasets`.
Also, you can download them from [here](https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/base_small.csr.gz) and [here](https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/queries.csr.gz).

For your application, you should also transform your dataset and queries into this binary format, then rename them as `[datasetName].csr` and `[queries].dev.csr`. Then, put them in the directory `./datasetss`.

### Download the datasets
There are 3 datasets that can be downloaded from the internet, which are `base_small`, `base_1M` and `base_full`. `base_1M` and `base_full` are used in the experiments in our paper (`SPLADE1M` and `SPLADE-Full`). To download them, we can follow the guidance below (Take `base_1M` as example):
#### Windows
We can download `base_1M` and the query sets directly from [here](https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/base_1M.csr.gz) and [here](https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/queries.csr.gz) via browser and then unzip them.

#### Linux
```bash
mkdir datasets
cd datasets
rm *
wget https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/queries.dev.csr.gz
wget https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/base_1M.csr.gz
gunzip *.gz
cd ..
```

## Reproduction in Linux
To quickly reproduce the results in our paper, a one-click shell `run.sh` is prepared. We can run it as following:
```bash
cd .
./run.sh datasetID
```
`datasetID` can be `0`, `1` or `2`, which represents the datasets `base_small`, `base_1M` and `base_full`, respectively.
## Results
The experimental result is saved in the directory `./results/` as the file `Running_result.txt`.

## Reference
[Efficient Approximate Maximum Inner Product Search over Sparse Vectors (Submitted to ICDE 2024)](https://github.com/Jacyhust/SOSIA/tree/main/Report)
