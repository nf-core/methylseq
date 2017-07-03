#!/usr/bin/env bash

data_path="/tmp"
if [ -d "./test_data" ]
then
    data_path="./test_data"
    echo "Found data directory in current working directory, using ./test_data/"
fi

curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
docker -v >/dev/null 2>&1 || { echo >&2 "I require docker, but it's not installed. Visit https://www.docker.com/products/overview#/install_the_platform  ."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

data_dir=${data_path}/ngi-bisulfite_test_set
if [ -d $data_dir ]
then
    echo "Found existing test set, using $data_dir"
else
    echo "Downloading test set..."
    curl https://export.uppmax.uu.se/b2013064/test-data/ngi-bisulfite_test_set.tar.bz2 > ${data_path}/ngi-bisulfite_test_set.tar.bz2
    echo "Unpacking test set..."
    tar xvjf ${data_path}/ngi-bisulfite_test_set.tar.bz2 -C ${data_path}
    echo "Done"
fi

# Do we build a reference genome index or just use a pre-existing one?
if [ -z $1 ]
then
    buildrefs="--saveReference --fasta ${data_dir}/references/WholeGenomeFasta/genome.fa"
else
    buildrefs="--bismark_index ${data_dir}/references/BismarkIndex/"
fi

# Detect Travis fork for dockerhub image if we can
if [ -z "$TRAVIS_REPO_SLUG" ]; then
    dockerfl=""
else
    dockerimg=$(echo "$TRAVIS_REPO_SLUG" | awk '{print tolower($0)}')
    echo "Detected repo as '$TRAVIS_REPO_SLUG' - using docker image '$dockerimg'"
    dockerfl="-with-docker $dockerimg"
fi

# Run name
run_name="Test MethylSeq Run: "$(date +%s)

cmd="nextflow run ../bismark.nf -resume -name \"$run_name\" -profile testing $dockerfl $buildrefs --singleEnd --reads \"${data_dir}/*.fastq.gz\" & sleep 540 ; kill \$!"
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
echo "-------------------------------------------------------"
echo "-------------------------------------------------------"
cat .nextflow.log

