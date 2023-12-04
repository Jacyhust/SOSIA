cd datasets
rm *

wget https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/queries.dev.csr.gz

wget https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/base_small.csr.gz
# wget https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/base_1M.csr.gz
# wget https://storage.googleapis.com/ann-challenge-sparse-vectors/csr/base_full.csr.gz
gunzip *.gz