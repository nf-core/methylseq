#!/usr/bin/env bash

# print_usage()
function print_usage {
  echo -e  "\nUsage:\t$0\n" \
    "\t\t[-a (aligner to use)\n" \
    "\t\t[-b (build genome references)\n" \
    "\t\t[-r (run in RRBS mode)\n" \
    "\t\t[-n (run in notrim mode)\n" \
    "\t\t[-p (run bwameth pipeline)\n" \
    "\t\t[-u (run UPPMAX test)\n" \
    "\t\t[-t <test data directory>]\n" \
    "\t\t[-d <docker image>]\n" \
    "\t\t[-s <singularity image>]\n" \
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
else
    dockerfl="-with-docker"
fi

# Look for an existing test data directory
data_path="./test_data"
data_dir=${data_path}/ngi-bisulfite_test_set

if [ -d $data_dir ]
then
    echo "Found existing test set, using $data_dir"
else
    echo "Unpacking test set..."
    tar xvjf ${data_path}/ngi-bisulfite_test_set.tar.bz2 -C ${data_path}
    echo "Done"
fi

# command line options
pipelinescript="../main.nf"
aligner="bismark"
profile="--max_cpus 2 --max_memory '7.GB' --max_time '48.h'"
if [ -d "${data_dir}/references/BismarkIndex/" ]
then
    refs="--bismark_index ${data_dir}/references/BismarkIndex/"
else
    refs="--bismark_index results/reference_genome/BismarkIndex/"
    echo "Couldn't find Bismark ref index in '${data_dir}/references/BismarkIndex/' - assuming one has been built by a test already"
fi
rrbs=""
notrim=""
singularityfl=""

while getopts ":brnpuht:d:s:a:" opt; do
  case $opt in
    a)
      echo "Using aligner $OPTARG" >&2
      aligner=$OPTARG
      ;;
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
      profile="-profile uppmax_devel"
      ;;
    t)
      echo "Test data path specified" >&2
      data_path=$OPTARG
      ;;
    d)
      echo "Using docker image $OPTARG" >&2
      dockerfl="-with-docker $OPTARG"
      ;;
    s)
      echo "Using singularity image $OPTARG" >&2
      singularityfl="-with-singularity $OPTARG"
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

cmd="nextflow run $pipelinescript -resume --aligner $aligner $profile $notrim $rrbs $dockerfl $singularityfl $refs --singleEnd --reads \"${data_dir}/SRR389222_sub*.fastq.gz\""
echo "Starting nextflow... Command:"
echo $cmd
echo "-------------------------------------------------------"
eval $cmd
