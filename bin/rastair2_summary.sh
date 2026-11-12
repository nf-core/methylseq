#!/usr/bin/env bash
#
# Summarise a rastair version 2 CpG BED (`rastair call --bed`) as MultiQC custom content.
#
# Usage: rastair2_summary.sh <sample> <bed.gz> [min_depth]
#
# Writes, for one sample:
#   <sample>.rastair2_generalstats_mqc.tsv  CpG methylation, CpGs covered and mean depth over chromosomes (General Statistics)
#   <sample>.rastair2_beta_hist_mqc.tsv     distribution of per-CpG methylation for CpGs with >= min_depth reads
#   <sample>.rastair2_contigs_mqc.tsv       CpG methylation per contig (chromosomes and any spike-ins or controls)
#
# Only reference CpGs (`cpg` column == REF) are counted. Methylation is pooled over reads:
# sum(mod) / sum(mod + unmod). Genome-wide figures use chromosome-like contigs (chr1, 2, chrX,
# chr1_KI270706v1_random, ...) so that spike-ins or controls do not skew them; if the reference has
# no chromosome-like contigs, all contigs are used. The per-contig table lists every contig with at
# least 10 covered CpGs, except unplaced, alternate, decoy, HLA and EBV sequences.

set -euo pipefail

sample=$1
bed=$2
min_depth=${3:-5}

gzip -cd "$bed" | awk -F '\t' -v sample="$sample" -v min_depth="$min_depth" '
    function is_minor(c) {
        return c ~ /(_alt|_random|_fix|_decoy|_hap[0-9]*)$|^chrUn|^HLA-|EBV/
    }
    function sort_key(c,    b) {
        b = c; sub(/^chr/, "", b)
        if (b ~ /^[0-9]+$/) return sprintf("0%06d", b + 0)
        if (b == "X") return "1X"
        if (b == "Y") return "2Y"
        if (b == "M" || b == "MT") return "3M"
        return "4" c
    }
    function is_chrom(c) {
        return c ~ /^(chr)?([0-9]+|[XYMW]|MT|Z)$/ || is_minor(c)
    }
    BEGIN {
        for (i = 0; i < 10; i++) { hist[i] = 0; ghist[i] = 0 }
    }
    /^#/ { next }
    $15 != "REF" { next }
    {
        unmod = $7; mod = $8; depth = unmod + mod
        if (depth == 0) next
        sites++; mod_sum += mod; depth_sum += depth
        contig_mod[$1] += mod; contig_depth[$1] += depth; contig_sites[$1]++
        bin = int(10 * mod / depth); if (bin > 9) bin = 9
        if (depth >= min_depth) hist[bin]++
        if (is_chrom($1)) {
            gsites++; gmod_sum += mod; gdepth_sum += depth
            if (depth >= min_depth) ghist[bin]++
        }
    }
    END {
        # Genome-wide figures over chromosomes only, unless there are none
        if (gsites > 0) {
            sites = gsites; mod_sum = gmod_sum; depth_sum = gdepth_sum
            for (i = 0; i < 10; i++) hist[i] = ghist[i]
        }
        gs = sample ".rastair2_generalstats_mqc.tsv"
        print "# id: \"rastair2_generalstats\"" > gs
        print "# plot_type: \"generalstats\"" > gs
        print "# headers:" > gs
        print "#     rastair2_cpg_meth:" > gs
        print "#         title: \"CpG meth\"" > gs
        print "#         namespace: \"rastair2\"" > gs
        print "#         description: \"CpG methylation pooled over reads, reference CpGs on chromosomes (rastair2)\"" > gs
        print "#         suffix: \"%\"" > gs
        print "#         min: 0" > gs
        print "#         max: 100" > gs
        print "#         scale: \"RdYlBu-rev\"" > gs
        print "#         format: \"{:,.1f}\"" > gs
        print "#     rastair2_cpg_sites:" > gs
        print "#         title: \"CpGs\"" > gs
        print "#         namespace: \"rastair2\"" > gs
        print "#         description: \"Reference CpG positions (per strand) with at least one informative read\"" > gs
        print "#         scale: \"Greens\"" > gs
        print "#         format: \"{:,.0f}\"" > gs
        print "#     rastair2_cpg_depth:" > gs
        print "#         title: \"CpG depth\"" > gs
        print "#         namespace: \"rastair2\"" > gs
        print "#         description: \"Mean informative reads per covered CpG position\"" > gs
        print "#         suffix: \"X\"" > gs
        print "#         scale: \"Blues\"" > gs
        print "#         format: \"{:,.1f}\"" > gs
        print "Sample\trastair2_cpg_meth\trastair2_cpg_sites\trastair2_cpg_depth" > gs
        printf "%s\t%.2f\t%d\t%.2f\n", sample, (depth_sum ? 100 * mod_sum / depth_sum : 0), sites, (sites ? depth_sum / sites : 0) > gs

        hf = sample ".rastair2_beta_hist_mqc.tsv"
        print "# id: \"rastair2_beta_hist\"" > hf
        print "# section_name: \"rastair2 CpG methylation distribution\"" > hf
        print "# description: \"Reference CpG positions on chromosomes in each methylation bin, for positions with at least " min_depth " informative reads.\"" > hf
        print "# plot_type: \"bargraph\"" > hf
        print "# pconfig:" > hf
        print "#     id: \"rastair2_beta_hist_plot\"" > hf
        print "#     title: \"rastair2: CpG methylation distribution\"" > hf
        print "#     ylab: \"CpG positions\"" > hf
        print "#     cpswitch_c_active: false" > hf
        header = "Sample"; line = sample
        for (i = 0; i < 10; i++) {
            header = header "\t" (10 * i) "-" (10 * i + 10) "%"
            line = line "\t" hist[i]
        }
        print header > hf
        print line > hf

        # Per-contig table; MultiQC rejects a table without columns, so skip it when no contig qualifies.
        # Columns in natural order: chr1..chr22, X, Y, M/MT, then other contigs (spike-ins) by name.
        n_contigs = 0
        for (c in contig_depth) {
            if (is_minor(c) || contig_sites[c] < 10) continue
            n_contigs++; keep[n_contigs] = c; key[n_contigs] = sort_key(c)
        }
        for (i = 2; i <= n_contigs; i++) {
            c = keep[i]; k = key[i]
            for (j = i - 1; j >= 1 && key[j] > k; j--) { keep[j + 1] = keep[j]; key[j + 1] = key[j] }
            keep[j + 1] = c; key[j + 1] = k
        }
        header = "Sample"; line = sample
        for (i = 1; i <= n_contigs; i++) {
            c = keep[i]
            header = header "\t" c
            line = line "\t" sprintf("%.2f", 100 * contig_mod[c] / contig_depth[c])
        }
        if (n_contigs > 0) {
            cf = sample ".rastair2_contigs_mqc.tsv"
            print "# id: \"rastair2_contigs\"" > cf
            print "# section_name: \"rastair2 CpG methylation per contig\"" > cf
            print "# description: \"CpG methylation (%) pooled over reads for each contig, including any methylated or unmethylated spike-in controls. Contigs with fewer than 10 covered CpG positions, and unplaced, alternate, decoy, HLA and EBV sequences, are left out.\"" > cf
            print "# plot_type: \"table\"" > cf
            print "# pconfig:" > cf
            print "#     id: \"rastair2_contigs_table\"" > cf
            print "#     title: \"rastair2: CpG methylation per contig\"" > cf
            print "#     min: 0" > cf
            print "#     max: 100" > cf
            print "#     format: \"{:,.1f}\"" > cf
            print header > cf
            print line > cf
        }
    }
'
