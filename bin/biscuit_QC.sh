#!/usr/bin/env bash
## make sure the following is in PATH
## biscuit samtools, bedtools, awk

# Use python's argparse module in shell scripts
#
# The function `argparse` parses its arguments using
# argparse.ArgumentParser; the parser is defined in the function's
# stdin.
#
# Executing ``argparse.bash`` (as opposed to sourcing it) prints a
# script template.
#
# https://github.com/nhoffman/argparse-bash
# MIT License - Copyright (c) 2015 Noah Hoffman

argparse(){
    argparser=$(mktemp 2>/dev/null || mktemp -t argparser)
    cat > "$argparser" <<EOF
from __future__ import print_function
import sys
import argparse
import os
class MyArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        """Print help and exit with error"""
        super(MyArgumentParser, self).print_help(file=file)
        sys.exit(1)
parser = MyArgumentParser(prog=os.path.basename("$0"),
            description="""$ARGPARSE_DESCRIPTION""")
EOF

    # stdin to this function should contain the parser definition
    cat >> "$argparser"

    cat >> "$argparser" <<EOF
args = parser.parse_args()
for arg in [a for a in dir(args) if not a.startswith('_')]:
    key = arg ####my change
    value = getattr(args, arg, None)
    if isinstance(value, bool) or value is None:
        print('{0}="{1}";'.format(key, 'yes' if value else ''))
    elif isinstance(value, list):
        print('{0}=({1});'.format(key, ' '.join('"{0}"'.format(s) for s in value)))
    else:
        print('{0}="{1}";'.format(key, value))
EOF

    # Define variables corresponding to the options if the args can be
    # parsed without errors; otherwise, print the text of the error
    # message.
    if python "$argparser" "$@" &> /dev/null; then
        eval $(python "$argparser" "$@")
        retval=0
    else
        python "$argparser" "$@"
        retval=1
    fi

    rm "$argparser"
    return $retval
}

#!/usr/bin/env bash
################################################################################
##
## Quality Control script for BISCUIT output
##
## Output from this script can be fed into MultiQC to produce a nice HTML output
## showing the different BISCUIT QC metrics
##
## Notes:
##   1.) biscuit, samtools, bedtools, and awk all must be in PATH for script to
##       work
##
## Created by:
##   Wanding Zhou
##
## Creation date:
##   May 2019
##
## Update notes:
##   Dec 2019 -
##     - Clean up code to make more readable
##     - Catch empty files, alert user, and remove files
##
################################################################################

# Check for biscuit, samtools, bedtools, awk in PATH
function check_path {
  if [[ `which biscuit 2>&1 > /dev/null` ]]; then
      echo "biscuit does not exist in PATH"
      exit 1
  fi
  if [[ `which samtools 2>&1 > /dev/null` ]]; then
      echo "samtools does not exist in PATH"
      exit 1
  fi
  if [[ `which bedtools 2>&1 > /dev/null` ]]; then
      echo "bedtools does not exist in PATH"
      exit 1
  fi
  if [[ `which awk 2>&1 > /dev/null` ]]; then
      echo "awk does not exist in PATH"
      exit 1
  fi
  if [[ `which python 2>&1 > /dev/null` ]]; then
      echo "python does not exist in PATH"
      exit 1
  fi
}

# Check that certain variables have been set and files exist
#TODO: Change "<unset>" to "NULL"/"NA"/something similar
#TODO: Also change in BISCUIT QC setup files
function check_variables {
    VARS="
    BISCUIT_CPGBED
    BISCUIT_CGIBED
    BISCUIT_RMSK
    BISCUIT_EXON
    BISCUIT_GENE
    BISCUIT_TOPGC_BED
    BISCUIT_BOTGC_BED
    input_bam
    input_vcf
    "

    for var in $VARS; do
        if [[ ${!var} != "<unset>" && ! -f ${!var} ]]; then
            >&2 echo "$var: ${!var} does not exist"
            exit 1
        fi
    done
}

# Check if QC files have at least some information
function basic_check_output_filled {
    prepend_path=$QCdir/${sname}
    TWO_LINE_FILES="
    ${prepend_path}_all_cv_table.txt
    ${prepend_path}_covdist_cpg_q40_botgc_table.txt
    ${prepend_path}_covdist_cpg_q40_table.txt
    ${prepend_path}_covdist_cpg_q40_topgc_table.txt
    ${prepend_path}_covdist_cpg_table.txt
    ${prepend_path}_covdist_q40_botgc_table.txt
    ${prepend_path}_covdist_q40_table.txt
    ${prepend_path}_covdist_q40_topgc_table.txt
    ${prepend_path}_covdist_table.txt
    ${prepend_path}_cpg_cv_table.txt
    ${prepend_path}_cpg_dist_table.txt
    ${prepend_path}_CpGRetentionByReadPos.txt
    ${prepend_path}_CpGRetentionDist.txt
    ${prepend_path}_CpHRetentionByReadPos.txt
    ${prepend_path}_freqOfTotalRetentionPerRead.txt
    ${prepend_path}_isize_score_table.txt
    ${prepend_path}_mapq_table.txt
    ${prepend_path}_totalBaseConversionRate.txt
    ${prepend_path}_totalReadConversionRate.txt
    "
    ONE_LINE_FILES="
    ${prepend_path}_dup_report.txt
    ${prepend_path}_strand_table.txt
    "

    echo "Running basic check on BISCUIT QC output"
    echo "Will remove any files that were obviously not filled properly"
    echo "This avoids clashes when running MultiQC"

    # All files that have a description line, then a table header line
    for FILE in ${TWO_LINE_FILES}; do
        if [[ ! -f "${FILE}" ]]; then
            >&2 echo "--- {FILE} --- was not initially created. Skipping!"
            continue
        fi
        if [[ `wc -l ${FILE} | awk '{print $1}'` -lt 3 ]]; then
            >&2 echo "--- ${FILE} --- has no entries. Check related files!"
            >&2 echo "Deleting --- ${FILE} --- since there are no entries to help with debugging."
            rm -f ${FILE}
        fi
    done

    # Files with only a description line
    for FILE in ${ONE_LINE_FILES}; do
        if [[ ! -f "${FILE}" ]]; then
            >&2 echo "--- {FILE} --- was not initially created. Skipping!"
            continue
        fi
        if [[ `wc -l ${FILE} | awk '{print $1}'` -lt 2 ]]; then
            >&2 echo "--- ${FILE} --- has no entries. Check related files!"
            >&2 echo "Deleting --- ${FILE} --- since there are no entries to help with debugging."
            rm -f ${FILE}
        fi
    done
}

function biscuitQC {

   # Simple check for necessary command line tools
    check_path

    # Check variables and their associated files exist
    check_variables

    # Create $QCdir if it does not exist
    if [ ! -d $QCdir ]; then
        mkdir -p $QCdir
    fi

    echo "Running BISCUIT QC"
    set -xe pipefail
  ##########################
  ## base coverage
  ##########################
  if [[ "$BISCUIT_QC_BASECOV" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_BASECOV ----"
   # bedtools genomecov -bga -split -ibam $input_bam -g ${BISCUIT_REFERENCE}.fai | bedtools sort >$QCdir/${sname}_bga.bed
   # samtools view -q 40 -b $input_bam | bedtools genomecov -ibam stdin -g ${BISCUIT_REFERENCE}.fai -bga -split | bedtools sort  >$QCdir/${sname}_bga_q40.bed
    bedtools genomecov -bga -split -ibam $input_bam -g ${BISCUIT_REFERENCE}.fai | LC_ALL=C sort --parallel=$processes  -k1,1 -k2,2n -T $QCdir >$QCdir/${sname}_bga.bed
    samtools view -q 40 -b $input_bam | bedtools genomecov -ibam stdin -g ${BISCUIT_REFERENCE}.fai -bga -split | LC_ALL=C sort --parallel=$processes   -k1,1 -k2,2n -T $QCdir  >$QCdir/${sname}_bga_q40.bed

    echo -e "BISCUITqc Depth Distribution (All)" >$QCdir/${sname}_covdist_table.txt
    echo -e "depth\tcount" >>$QCdir/${sname}_covdist_table.txt
    awk '{cnt[$4]+=$3-$2}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_bga.bed | sort -k1,1n -T $QCdir >>$QCdir/${sname}_covdist_table.txt

    echo -e "BISCUITqc Depth Distribution (Q40)" >$QCdir/${sname}_covdist_q40_table.txt
    echo -e "depth\tcount" >>$QCdir/${sname}_covdist_q40_table.txt
    awk '{cnt[$4]+=$3-$2}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_bga_q40.bed | sort -k1,1n -T $QCdir >>$QCdir/${sname}_covdist_q40_table.txt
  fi
  
  ##########################
  ## duplicate_coverage
  ##########################
  [[ ! -f "$QCdir/${sname}_bga.bed" ]] && BISCUIT_QC_DUPLICATE=false
  [[ ! -f "$QCdir/${sname}_bga_q40.bed" ]] && BISCUIT_QC_DUPLCIATE=false
  if [[ "$BISCUIT_QC_DUPLICATE" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_DUPLICATE ----"
    # duplicate
    #samtools view -f 0x400 -b $input_bam | bedtools genomecov -ibam stdin -g $BISCUIT_REFERENCE.fai -bga -split | bedtools sort >$QCdir/${sname}_bga_dup.bed
    samtools view -f 0x400 -b $input_bam | bedtools genomecov -ibam stdin -g $BISCUIT_REFERENCE.fai -bga -split | LC_ALL=C sort --parallel=$processes -k1,1 -k2,2n -T $QCdir >$QCdir/${sname}_bga_dup.bed

    # duplication rate
    echo -e "BISCUITqc Read Duplication Table" >$QCdir/${sname}_dup_report.txt
    echo -ne "#bases covered by all reads: " >>$QCdir/${sname}_dup_report.txt
    awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga.bed >>$QCdir/${sname}_dup_report.txt
    echo -ne "#bases covered by duplicate reads: " >>$QCdir/${sname}_dup_report.txt
    awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga_dup.bed >>$QCdir/${sname}_dup_report.txt

    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      # high GC content
      echo -ne "#high-GC bases covered by all reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga.bed -b $BISCUIT_TOPGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#high-GC bases covered by duplicate reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup.bed -b $BISCUIT_TOPGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt

      # low GC content
      echo -ne "#low-GC bases covered by all reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga.bed -b $BISCUIT_BOTGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#low-GC bases covered by duplicate reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup.bed -b $BISCUIT_BOTGC_BED -sorted | awk 'BEGIN{a=0}$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
    fi
    
    ## Q40
    # duplicate
   # samtools view -f 0x400 -q 40 -b $input_bam | bedtools genomecov -ibam stdin -g $BISCUIT_REFERENCE.fai -bga -split | bedtools sort >$QCdir/${sname}_bga_dup_q40.bed
    samtools view -f 0x400 -q 40 -b $input_bam | bedtools genomecov -ibam stdin -g $BISCUIT_REFERENCE.fai -bga -split | LC_ALL=C sort --parallel=$processes -k1,1 -k2,2n -T $QCdir >$QCdir/${sname}_bga_dup_q40.bed

    # duplication rate
    echo -ne "#bases covered by all q40-reads: " >>$QCdir/${sname}_dup_report.txt
    awk '$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga_q40.bed >>$QCdir/${sname}_dup_report.txt
    echo -ne "#bases covered by duplicate q40-reads: " >>$QCdir/${sname}_dup_report.txt
    awk '$4>0{a+=$3-$2}END{print a}' $QCdir/${sname}_bga_dup_q40.bed >>$QCdir/${sname}_dup_report.txt

    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      # high GC content
      echo -ne "#high-GC bases covered by all q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#high-GC bases covered by duplicate q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt

      # low GC content
      echo -ne "#low-GC bases covered by all q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
      echo -ne "#low-GC bases covered by duplicate q40-reads: " >>$QCdir/${sname}_dup_report.txt
      bedtools intersect -a $QCdir/${sname}_bga_dup_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk '$4>0{a+=$3-$2}END{print a}' >>$QCdir/${sname}_dup_report.txt
    fi
  fi
  
  ##########################
  ## cpg coverage
  ##########################

  [[ ! -f "$BISCUIT_CPGBED" ]] && BISCUIT_QC_CPGCOV=false
  [[ ! -f "$QCdir/${sname}_bga.bed" ]] && BISCUIT_QC_CPGCOV=false
  [[ ! -f "$QCdir/${sname}_bga_q40.bed" ]] && BISCUIT_QC_CPGCOV=false
  if [[ "$BISCUIT_QC_CPGCOV" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_CPGCOV ----"
    bedtools intersect -a $BISCUIT_CPGBED -b $QCdir/${sname}_bga.bed -wo -sorted | bedtools groupby -g 1-3 -c 7 -o min >$QCdir/${sname}_cpg.bed
    bedtools intersect -a $BISCUIT_CPGBED -b $QCdir/${sname}_bga_q40.bed -wo -sorted | bedtools groupby -g 1-3 -c 7 -o min >$QCdir/${sname}_cpg_q40.bed

    echo -e "BISCUITqc CpG Depth Distribution (All)" >$QCdir/${sname}_covdist_cpg_table.txt
    echo -e "depth\tcount" >>$QCdir/${sname}_covdist_cpg_table.txt
    awk '{cnt[$4]+=1}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_cpg.bed | sort -k1,1n >>$QCdir/${sname}_covdist_cpg_table.txt

    echo -e "BISCUITqc CpG Depth Distribution (Q40)" >$QCdir/${sname}_covdist_cpg_q40_table.txt
    echo -e "depth\tcount" >>$QCdir/${sname}_covdist_cpg_q40_table.txt
    awk '{cnt[$4]+=1}END{for(cov in cnt) {print int(cov)"\t"int(cnt[cov]);}}' $QCdir/${sname}_cpg_q40.bed | sort -k1,1n >>$QCdir/${sname}_covdist_cpg_q40_table.txt
  fi

  ##########################
  ## cpg distribution
  ##########################

  [[ ! -f "$QCdir/${sname}_cpg_q40.bed" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$QCdir/${sname}_cpg.bed" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_EXON" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_RMSK" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_GENE" ]] && BISCUIT_QC_CPGDIST=false
  [[ ! -f "$BISCUIT_CGIBED" ]] && BISCUIT_QC_CPGDIST=false
  if [[ "$BISCUIT_QC_CPGDIST" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_CPGDIST ----"
    # whole genome
    echo -e "BISCUITqc CpG Distribution Table" >$QCdir/${sname}_cpg_dist_table.txt
    wc -l $QCdir/${sname}_cpg_q40.bed | awk -F" " '{printf("Territory\tAll\tUniqCov\tAllCov\nTotalCpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    awk '$4>0{a+=1}END{printf("\t%d",a)}' $QCdir/${sname}_cpg_q40.bed >>$QCdir/${sname}_cpg_dist_table.txt
    awk '$4>0{a+=1}END{printf("\t%d\n",a)}' $QCdir/${sname}_cpg.bed >>$QCdir/${sname}_cpg_dist_table.txt

    # exon
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_EXON) -sorted | wc -l | awk -F" " '{printf("ExonicCpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_EXON) -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b <(bedtools merge -i $BISCUIT_EXON) -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt

    # repeat
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_RMSK) -sorted | wc -l | awk -F" " '{printf("RepeatCpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_RMSK) -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b <(bedtools merge -i $BISCUIT_RMSK) -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt

    # gene
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_GENE) -sorted | wc -l | awk -F" " '{printf("GenicCpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_GENE) -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b <(bedtools merge -i $BISCUIT_GENE) -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt

    # CGI
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_CGIBED) -sorted | wc -l | awk -F" " '{printf("CGICpGs\t%s",$1)}' >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_CGIBED) -sorted | awk '$4>0{a+=1}END{printf("\t%d",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg.bed -b <(bedtools merge -i $BISCUIT_CGIBED) -sorted | awk '$4>0{a+=1}END{printf("\t%d\n",a)}'  >>$QCdir/${sname}_cpg_dist_table.txt

    >&2 echo "`date`---- BISCUIT_QC_CGICOV ----"
    # how CGI is covered by at least one q40-read in at least one CpG
    echo >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands\t" >>$QCdir/${sname}_cpg_dist_table.txt
    zcat $BISCUIT_CGIBED | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b <(bedtools merge -i $BISCUIT_CGIBED) -sorted -wo | awk '$4>0{print $5":"$6"-"$7}' | uniq -c | awk -F" " '{print $2"\t"$1}' >> $QCdir/${sname}_cpg_dist_table_temp.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least one CpG\t" >>$QCdir/${sname}_cpg_dist_table.txt
    less $QCdir/${sname}_cpg_dist_table_temp.txt | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least three CpGs\t" >>$QCdir/${sname}_cpg_dist_table.txt
    awk -F" " '$2>=3' $QCdir/${sname}_cpg_dist_table_temp.txt | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least five CpGs\t" >>$QCdir/${sname}_cpg_dist_table.txt
    awk -F" " '$2>=5' $QCdir/${sname}_cpg_dist_table_temp.txt | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    echo -ne "#CpG Islands covered by at least one q40-read in at least ten CpGs\t" >>$QCdir/${sname}_cpg_dist_table.txt
    awk -F" " '$2>=10' $QCdir/${sname}_cpg_dist_table_temp.txt | wc -l >>$QCdir/${sname}_cpg_dist_table.txt
    /bin/rm $QCdir/${sname}_cpg_dist_table_temp.txt
  fi
  
  ##########################
  ## uniformity
  ##########################
  [[ ! -f "$QCdir/${sname}_covdist_q40_table.txt" ]] && BISCUIT_QC_UNIFORMITY=false
  [[ ! -f "$QCdir/${sname}_bga_q40.bed" ]] && BISCUIT_QC_UNIFORMITY=false
  if [[ "$BISCUIT_QC_UNIFORMITY" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_UNIFORMITY ----"

    echo -e "BISCUITqc Uniformity Table" >$QCdir/${sname}_all_cv_table.txt
    awk -v sname="${sname}" '{cnt[$1]=$2}END{for (cov in cnt) {sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print "sample\tmu\tsigma\tcv\n"sname"_all\t"mu"\t"sigma"\t"sigma/mu}' $QCdir/${sname}_covdist_q40_table.txt >>$QCdir/${sname}_all_cv_table.txt

    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      echo -e "BISCUITqc Depth Distribution (high GC, Q40)" >$QCdir/${sname}_covdist_q40_topgc_table.txt
      echo -e "depth\tcount" >>$QCdir/${sname}_covdist_q40_topgc_table.txt
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk -v sname="${sname}" -v output="$QCdir/${sname}_all_cv_table.txt" '{cnt[$4]+=$3-$2}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_all_topgc\t"mu"\t"sigma"\t"sigma/mu >>output}' | sort -k1,1n >>$QCdir/${sname}_covdist_q40_topgc_table.txt

      echo -e "BISCUITqc Depth Distribution (low GC, Q40)" >$QCdir/${sname}_covdist_q40_botgc_table.txt
      echo -e "depth\tcount" >>$QCdir/${sname}_covdist_q40_botgc_table.txt
      bedtools intersect -a $QCdir/${sname}_bga_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk -v sname="${sname}" -v output="$QCdir/${sname}_all_cv_table.txt" '{cnt[$4]+=$3-$2}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_all_botgc\t"mu"\t"sigma"\t"sigma/mu >>output}' | sort -k1,1n >>$QCdir/${sname}_covdist_q40_botgc_table.txt
    fi
  fi
  
  ##########################
  ## cpg uniformity
  ##########################
  [[ ! -f "$QCdir/${sname}_covdist_cpg_q40_table.txt" ]] && BISCUIT_QC_CPGUNIF=false
  [[ ! -f "$QCdir/${sname}_cpg_q40.bed" ]] && BISCUIT_QC_CPGUNIF=false
  if [[ "$BISCUIT_QC_CPGUNIF" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_CPGUNIF ----"

    echo -e "BISCUITqc CpG Uniformity Table" >$QCdir/${sname}_cpg_cv_table.txt
    awk -v sname="${sname}" '{cnt[$1]=$2}END{for(cov in cnt) {sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print "sample\tmu\tsigma\tcv\n"sname"_cpg\t"mu"\t"sigma"\t"sigma/mu}' $QCdir/${sname}_covdist_cpg_q40_table.txt >>$QCdir/${sname}_cpg_cv_table.txt

    if [[ -f "$BISCUIT_TOPGC_BED" && -f "$BISCUIT_BOTGC_BED" ]]; then
      echo -e "BISCUITqc CpG Depth Distribution (high GC, Q40)" >$QCdir/${sname}_covdist_cpg_q40_topgc_table.txt
      echo -e "depth\tcount" >>$QCdir/${sname}_covdist_cpg_q40_topgc_table.txt
      bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_TOPGC_BED -sorted | awk -v sname="${sname}" -v output="$QCdir/${sname}_cpg_cv_table.txt" '{cnt[$4]+=1}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_cpg_topgc\t"mu"\t"sigma"\t"sigma/mu >>output}' | sort -k1,1n >>$QCdir/${sname}_covdist_cpg_q40_topgc_table.txt

      echo -e "BISCUITqc CpG Depth Distribution (low GC, Q40)" >$QCdir/${sname}_covdist_cpg_q40_botgc_table.txt
      echo -e "depth\tcount" >>$QCdir/${sname}_covdist_cpg_q40_botgc_table.txt
      bedtools intersect -a $QCdir/${sname}_cpg_q40.bed -b $BISCUIT_BOTGC_BED -sorted | awk -v sname="${sname}" -v output="$QCdir/${sname}_cpg_cv_table.txt" '{cnt[$4]+=1}END{for (cov in cnt) {print cov"\t"cnt[cov]; sum_cov+=cnt[cov]*cov; sum_cnt+=cnt[cov];} for(cov in cnt) {sum_var+=((cov-mu)^2)*cnt[cov];} mu=sum_cov/sum_cnt; sigma=sqrt(sum_var/sum_cnt); print sname"_cpg_botgc\t"mu"\t"sigma"\t"sigma/mu >>output}' | sort -k1,1n >>$QCdir/${sname}_covdist_cpg_q40_botgc_table.txt
    fi
  fi
  
  ##########################
  ## bisulfite conversion
  ##########################
  [[ ! -f "$input_vcf" ]] && BISCUIT_QC_BSCONV=false
  if [[ "$BISCUIT_QC_BSCONV" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_BSCONV ----"
    
	#echo -e "BISCUITqc Frequency of Total Retention per Read Table" >$QCdir/${sname}_freqOfTotalRetentionPerRead.txt
    #samtools view -h -q 40 $input_bam | biscuit bsconv $BISCUIT_REFERENCE - | awk 'match($0,/ZN:Z:([^ ]*)/,a){print gensub(/[A-Z,_]+/, "\t", "g", a[1])}' | cut -f2,4,6,8 | awk -v OFS="\t" '{ra[$1]+=1;rc[$2]+=1;rg[$3]+=1;rt[$4]+=1;}END{for(k in ra) {print "CA", k, ra[k]} for(k in rc) {print "CC", k, rc[k]} for(k in rg) {print "CG", k, rg[k]} for(k in rt) {print "CT", k, rt[k]}}' | sort -k1,1 -k2,2n | awk 'BEGIN{print "CTXT\tnumRET\tCnt"}{print}' >>$QCdir/${sname}_freqOfTotalRetentionPerRead.txt
	
	echo -e "BISCUITqc Frequency of Total Retention per Read Table\nCTXT\tnumRET\tCnt" >$QCdir/${sname}_freqOfTotalRetentionPerRead.txt
	samtools view -h -q 40 $input_bam | biscuit bsconv $BISCUIT_REFERENCE - |  grep -E -o --color "ZN:Z:[^ ].*" | awk -F '[^0-9]*' -v OFS="\t" '{ra[$2]+=1;rc[$4]+=1;rg[$6]+=1;rt[$8]+=1;}END{for(k in ra) {print "CA", k, ra[k]} for(k in rc) {print "CC", k, rc[k]} for(k in rg) {print "CG", k, rg[k]} for(k in rt) {print "CT", k, rt[k]}}' | sort -k1,1 -k2,2n  >>$QCdir/${sname}_freqOfTotalRetentionPerRead.txt
    
	echo -e "BISCUITqc Conversion Rate by Base Average Table" >$QCdir/${sname}_totalBaseConversionRate.txt
    biscuit vcf2bed -et c $input_vcf | awk '{beta_sum[$6]+=$8; beta_cnt[$6]+=1;} END{print "CA\tCC\tCG\tCT"; print beta_sum["CA"]/beta_cnt["CA"]"\t"beta_sum["CC"]/beta_cnt["CC"]"\t"beta_sum["CG"]/beta_cnt["CG"]"\t"beta_sum["CT"]/beta_cnt["CT"];}' >>$QCdir/${sname}_totalBaseConversionRate.txt

    echo -e "BISCUITqc Conversion Rate by Read Average Table" >$QCdir/${sname}_totalReadConversionRate.txt
    samtools view -hq 40 -F 0x900 $input_bam | biscuit bsconv -b $BISCUIT_REFERENCE - | awk '{for(i=1;i<=8;++i) a[i]+=$i;}END{print "CpA\tCpC\tCpG\tCpT"; print a[1]/(a[1]+a[2])"\t"a[3]/(a[3]+a[4])"\t"a[5]/(a[5]+a[6])"\t"a[7]/(a[7]+a[8]);}' >>$QCdir/${sname}_totalReadConversionRate.txt

    echo -e "BISCUITqc CpH Retention by Read Position Table" >$QCdir/${sname}_CpHRetentionByReadPos.txt
    echo -e "ReadInPair\tPosition\tConversion/Retention\tCount" >>$QCdir/${sname}_CpHRetentionByReadPos.txt
    #samtools view -hq 40 $input_bam | biscuit cinread $BISCUIT_REFERENCE - -t ch -p QPAIR,CQPOS,CRETENTION | sort | uniq -c | awk -F" " '$4!="N"{print $2"\t"$3"\t"$4"\t"$1}' | sort -k1,1 -k2,2n -T $QCdir >>$QCdir/${sname}_CpHRetentionByReadPos.txt
    samtools view -hq 40 $input_bam | biscuit cinread $BISCUIT_REFERENCE - -t ch -p QPAIR,CQPOS,CRETENTION | awk -v OFS="\t"  '$3!="N"{sorting[$0]++ } END {for (i in sorting) print i,sorting[i]}' |  sort -k1,1 -k2,2n >>$QCdir/${sname}_CpHRetentionByReadPos.txt

    echo -e "BISCUITqc CpG Retention by Read Position Table" >$QCdir/${sname}_CpGRetentionByReadPos.txt
    echo -e "ReadInPair\tPosition\tConversion/Retention\tCount" >>$QCdir/${sname}_CpGRetentionByReadPos.txt
    #samtools view -hq 40 $input_bam | biscuit cinread $BISCUIT_REFERENCE - -t cg -p QPAIR,CQPOS,CRETENTION | sort | uniq -c | awk -F" " '$4!="N"{print $2"\t"$3"\t"$4"\t"$1}' | sort -k1,1 -k2,2n -T $QCdir >>$QCdir/${sname}_CpGRetentionByReadPos.txt
    samtools view -hq 40 $input_bam | biscuit cinread $BISCUIT_REFERENCE - -t cg -p QPAIR,CQPOS,CRETENTION | awk -v OFS="\t"  '$3!="N"{sorting[$0]++ } END {for (i in sorting) print i,sorting[i]}' |  sort -k1,1 -k2,2n >>$QCdir/${sname}_CpGRetentionByReadPos.txt
  fi

  ####################
  ## mapping_summary
  ####################
  if [[ "$BISCUIT_QC_MAPPING" == true ]]; then
    >&2 echo "`date`---- BISCUIT_QC_MAPPING ----"
    echo -e "BISCUITqc Strand Table" >$QCdir/${sname}_strand_table.txt
    biscuit cinread -p QPAIR,STRAND,BSSTRAND $BISCUIT_REFERENCE $input_bam | awk '{a[$1$2$3]+=1}END{for(strand in a) {print "strand\t"strand"\t"a[strand];}}' >>$QCdir/${sname}_strand_table.txt

    echo -e "BISCUITqc Mapping Quality Table" >$QCdir/${sname}_mapq_table.txt
    echo -e "MapQ\tCount" >>$QCdir/${sname}_mapq_table.txt
    samtools view -F 0x100 -f 0x4 $input_bam | wc -l | cat <(echo -ne "unmapped\t") - >>$QCdir/${sname}_mapq_table.txt
    samtools view -F 0x104 $input_bam | awk '{cnt[$5]+=1}END{for(mapq in cnt) {print mapq"\t"cnt[mapq];}}' | sort -k1,1n >>$QCdir/${sname}_mapq_table.txt
    ## insert size
    ## this excludes read by AS (40) and mapq (40)
    echo -e "BISCUITqc Insert Size, Score Table" >$QCdir/${sname}_isize_score_table.txt
    echo -e "InsertSize/Score\tValue\tFraction" >>$QCdir/${sname}_isize_score_table.txt
   # samtools view -F 0x104 $input_bam | awk '{match($0,/AS:i:([0-9]*)/,a); score[a[1]]+=1; sumscore+=1; if (and($2,0x2) && a[1]>=40 && $5>=40 && $9>=0 && $9 <=2000) {isize[$9]+=1; sumisize+=1}}END{for(k in isize){print "I\t"k"\t"isize[k] / sumisize} for(k in score){print "S\t"k"\t"score[k] / sumscore}}' | sort -k1,1 -k2,2n >>$QCdir/${sname}_isize_score_table.txt
    samtools view -F 0x104 $input_bam | grep -E 'AS:i:([0-9]*)' | awk -F 'AS:i:|\t' '{score[$17]+=1; sumscore+=1; if (and($2,0x2) && $17>=40 && $5>=40 && $9>=0 && $9 <=2000) {isize[$9]+=1; sumisize+=1}}END{for(k in isize){print "I\t"k"\t"isize[k] / sumisize} for(k in score){print "S\t"k"\t"score[k] / sumscore}}' | sort -k1,1 -k2,2n >>$QCdir/${sname}_isize_score_table.txt
  fi


  ###################################
  ## CpG retention distribution
  ###################################
  [[ ! -f "$input_vcf" ]] && BISCUIT_QC_BETAS=false
  if [[ "$BISCUIT_QC_BETAS" == true ]]; then
    echo -e "BISCUITqc Retention Distribution Table" >$QCdir/${sname}_CpGRetentionDist.txt
    echo -e "RetentionFraction\tCount" >>$QCdir/${sname}_CpGRetentionDist.txt
    biscuit vcf2bed -t cg $input_vcf | awk '$5>=3{a[sprintf("%3.0f", $4*100)]+=1}END{for (beta in a) print beta"\t"a[beta];}' | sort -k1,1n >>$QCdir/${sname}_CpGRetentionDist.txt
  fi
  
  ###################################
  ## Remove bed files 
  ###################################
  if [[ "$remove_bed_files" ]]; then
    rm $QCdir/*.bed
  fi
  
    ########################################
    ## Running check on output files
    ########################################
    basic_check_output_filled 
}




ARGPARSE_DESCRIPTION="Run QC on biscuit output"
argparse "$@" <<EOF || exit 1
parser.add_argument('assets_directory',help='path to assets directory file')
parser.add_argument('genome',help='path to genome.fa file')
parser.add_argument('sample_name',help='prefix of output files in the output directory')
parser.add_argument('input_bam',help='input bam file from alignment')
parser.add_argument('-v', '--vcf', type=str,default='<unset>',
                    help='Path to vcf outupt from BISCUIT')
parser.add_argument('-o', '--outdir', type=str, default='BISCUITqc',
                    help='output directory [default %(default)s]')
parser.add_argument('--do_not_remove_bed', action='store_false',
                    default=True, help='Whether remove bed files [default %(default)s]')
parser.add_argument('-p', '--processes', type=int, default=1,
                    help='Number of processes to use [default %(default)s]')
EOF

if [[ ! -d "$assets_directory" ]]; then
  echo "Assets directory missing: $assets_directory.";
  exit 1;
fi

source $(dirname ${BASH_SOURCE[0]})/setup.sh $genome $assets_directory 
input_vcf=$vcf
QCdir=$outdir
sname=$sample_name
remove_bed_files=$do_not_remove_bed

>&2 echo "## Running BISCUIT QC script with following configuration ##"
>&2 echo "=============="
>&2 echo "sample name: $sname"
>&2 echo "input bam:   $input_bam"
>&2 echo "input vcf:   $input_vcf"
>&2 echo "output dir:  $QCdir"
>&2 echo "REFERENCE:   $BISCUIT_REFERENCE"
>&2 echo "CPGBED:      $BISCUIT_CPGBED"
>&2 echo "CGIBED:      $BISCUIT_CGIBED"
>&2 echo "RMSK:        $BISCUIT_RMSK"
>&2 echo "EXON:        $BISCUIT_EXON"
>&2 echo "GENE:        $BISCUIT_GENE"
>&2 echo "TOPGC_BED:   $BISCUIT_TOPGC_BED"
>&2 echo "BOTGC_BED:   $BISCUIT_BOTGC_BED"
>&2 echo "=============="
biscuitQC
#>&2 echo $remove_bed_files $QCdir $input_vcf

>&2 echo -e "\nDone."
