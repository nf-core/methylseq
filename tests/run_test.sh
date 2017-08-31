#!/usr/bin/env bash

# print_usage()
function print_usage {
  echo -e  "\nUsage:\t$0\n" \
    "\t\t[-b (build genome references)\n" \
    "\t\t[-r (run in RRBS mode)\n" \
    "\t\t[-n (run in notrim mode)\n" \
    "\t\t[-p (run bwameth pipeline)\n" \
    "\t\t[-u (run UPPMAX test)\n" \
    "\t\t[-t <test data directory>]\n" \
    "\t\t[-d <docker image>]\n" \
    "\t\t[-h (show this help message)]" >&2 ;
}

# Check that we have required commands
curl --version >/dev/null 2>&1 || { echo >&2 "I require curl, but it's not installed. Aborting."; exit 1; }
tar --version >/dev/null 2>&1 || { echo >&2 "I require tar, but it's not installed. Aborting."; exit 1; }
nextflow -v >/dev/null 2>&1 || { echo >&2 "I require nextflow, but it's not installed. If you hava Java, run 'curl -fsSL get.nextflow.io | bash'. If not, install Java."; exit 1; }

# Detect Travis fork for dockerhub image if we can
dockerfl=""
if [[ ! -z "$TRAVIS_REPO_SLUG" ]]; then
    dockerimg=$(echo "$TRAVIS_REPO_SLUG" | awk '{print tolower($0)}')
    echo "Detected repo as '$TRAVIS_REPO_SLUG' - using docker image '$dockerimg'"
    dockerfl="-with-docker $dockerimg"
fi

# Look for an existing test data directory
data_path="/tmp"
if [ -d "./test_data" ]
then
    data_path="./test_data"
    echo "Found data directory in current working directory, using ./test_data/"
fi
data_dir=${data_path}/ngi-bisulfite_test_set

# command line options
pipelinescript="../bismark.nf"
profile="-profile testing"
refs="--bismark_index ${data_dir}/references/BismarkIndex/"
rrbs=""
notrim=""

while getopts ":brnpuht:d:" opt; do
  case $opt in
    b)
      echo "Building genome references" >&2
      refs="--saveReference --fasta ${data_dir}/references/WholeGenomeFasta/genome.fa"
      buildrefs=1
      ;;
    r)
      echo "Running in RRBS mode" >&2
      rrbs="--rrbs"
      ;;
    n)
      echo "Running in no-trimming mode" >&2
      notrim="--notrim"
      ;;
    p)
      echo "Running BWAmeth pipeline" >&2
      pipelinescript="../bwa-meth.nf --fasta ${data_dir}/references/WholeGenomeFasta/genome.fa"
      bwameth=1
      ;;
    u)
      echo "Running UPPMAX config" >&2
      profile="-profile devel"
      if [ data_path -ne "./test_data" ]; then
        data_path=$SNIC_NOBACKUP
      fi
      ;;
    t)
      echo "Test data path specified" >&2
      data_path=$OPTARG
      ;;
    d)
      echo "Using docker image $OPTARG" >&2
      dockerfl="-with-docker $OPTARG"
      ;;
    h)
      print_usage
      exit
      ;;
    :)
      echo -e "\nOption -$OPTARG requires an argument." >&2
      print_usage
      exit 1;
      ;;
    \?)
      print_usage
      echo -e "\nInvalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if [[ $buildrefs ]] && [[ $bwameth ]]; then
  refs="--saveReference --fasta_index ${data_dir}/references/WholeGenomeFasta/genome.fa.fai --bwa_meth_index results/reference_genome/genome.fa"
fi


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

# Run name
run_name="Test MethylSeq Run: "$(date +%s)

cmd="nextflow run $pipelinescript -resume -name \"$run_name\" $profile $notrim $rrbs $dockerfl $refs --singleEnd --reads \"${data_dir}/*.fastq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
