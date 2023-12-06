
# make clean
make

# echo $#
id=2
if [ $# -ge 1 ]
then
    id=$1
fi

datasets=( base_small base_1M base_full )
dpath=./datasets/${datasets[$id]}.csr
if [ -e ${dpath} ]
then
    echo The dataset ${datasets[$id]} is at ${dpath}
else
    echo The dataset ${datasets[$id]} does not exist. Downloading it now ...
    ./download.sh ${datasets[$id]}
fi

./sos 0 ${datasets[$id]} queries
./sos 1 ${datasets[$id]} queries
./sos 2 ${datasets[$id]} queries
./sos 3 ${datasets[$id]} queries
./sos 4 ${datasets[$id]} queries
