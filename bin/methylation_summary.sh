#!/usr/bin/env bash
#
# Summarise per-CpG methylation calls as MultiQC custom content.
#
# Usage: methylation_summary.sh [-f rastair2|bedgraph] [-n name] [-d min_depth] [-c control_contigs] [--] <sample> <file>...
#
#   -f  input format (default rastair2):
#         rastair2  one rastair version 2 CpG BED (`rastair call --bed`); columns found by name in the
#                   header (unmod, mod, cpg), and only reference CpGs (cpg == REF) counted
#         bedgraph  MethylDackel bedGraphs or Bismark coverage files (chrom, start, end, %, methylated,
#                   unmethylated), plain or gzipped. A file named *_CHG.* or *_CHH.* holds that context;
#                   any other holds CpGs.
#   -n  tool name for the MultiQC IDs and column names: rastair2 (default), methyldackel or bismark
#   -d  minimum informative reads for a CpG to count in the methylation distribution (default 5)
#   -c  comma-separated spike-in or control contigs (e.g. lambda,pUC19) to leave out of the genome-wide figures
#
# Writes, for one sample:
#   <sample>.<name>_generalstats_mqc.tsv  CpG methylation, CpGs covered and mean depth over the genome, plus
#                                         CHG and CHH methylation if given (General Statistics)
#   <sample>.<name>_beta_hist_mqc.tsv     distribution of per-CpG methylation for CpGs with >= min_depth reads
#   <sample>.<name>_contigs_mqc.tsv       CpG methylation per contig (chromosomes and any spike-ins or controls)
#
# Methylation is pooled over reads: sum(mod) / sum(mod + unmod). The genome-wide figures (General
# Statistics and the distribution) use chromosome-like contigs (chr1, 2, chrX, 2L, chrIV,
# chr1_KI270706v1_random, ...) if the input has any, otherwise every contig. Either way they leave out
# organelles (chrM, MT, chloroplast), decoys, EBV and the -c contigs, so that spike-ins and controls do
# not skew them. The contigs used are logged to stderr; if none of them has covered CpGs, the
# genome-wide files are not written. The per-contig table lists every contig with at least 10 covered
# CpGs (the 50 with the most, if there are more), except unplaced, alternate, decoy, HLA and EBV sequences.

set -euo pipefail
# Locale-independent number formatting (some awks write 60,37 under a comma-decimal LC_NUMERIC)
export LC_ALL=C

usage() {
    echo "Usage: methylation_summary.sh [-f rastair2|bedgraph] [-n name] [-d min_depth] [-c control_contigs] [--] <sample> <file>..." >&2
    exit 2
}

format=rastair2
name=rastair2
min_depth=5
# Passed through the environment: awk -v would interpret backslash escapes in the names
export control_contigs=""
while getopts "f:n:d:c:" opt; do
    case "$opt" in
        f) format=$OPTARG ;;
        n) name=$OPTARG ;;
        d) min_depth=$OPTARG ;;
        c) control_contigs=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))
[ $# -ge 2 ] || usage
[[ $format == rastair2 || $format == bedgraph ]] || { echo "methylation_summary.sh: -f must be rastair2 or bedgraph, got '$format'" >&2; exit 2; }
[[ $format == bedgraph || $# -eq 2 ]] || { echo "methylation_summary.sh: -f rastair2 takes one BED" >&2; exit 2; }
[[ $name =~ ^[a-z0-9_]+$ ]] || { echo "methylation_summary.sh: -n must be lower-case letters, digits and underscores, got '$name'" >&2; exit 2; }
[[ $min_depth =~ ^[0-9]+$ ]] || { echo "methylation_summary.sh: -d must be a whole number, got '$min_depth'" >&2; exit 2; }
case "$name" in
    methyldackel) label=MethylDackel ;;
    bismark) label=Bismark ;;
    *) label=$name ;;
esac
export sample=$1 name label
shift

# The bedgraph format has no context column: mark where each file starts with a "#context" line
for file in "$@"; do
    if [ "$format" = bedgraph ]; then
        case "$file" in
            *_CHG.*) printf '#context\tCHG\n' ;;
            *_CHH.*) printf '#context\tCHH\n' ;;
            *) printf '#context\tCpG\n' ;;
        esac
    fi
    gzip -cdf -- "$file"
done | awk -F '\t' -v min_depth="$min_depth" -v format="$format" '
    # Organelles: mitochondrion (chrM, MT, Mito, MtDNA, NC_012920.1) and chloroplast (chrC, Pt)
    function is_organelle(c,    lc) {
        lc = tolower(c)
        return lc ~ /^(chr)?(m|mt|mito|mtdna|mitochondrion(_genome)?|c|pt|chloroplast)$/ || c ~ /^NC_012920(\.[0-9]+)?$/
    }
    # Not part of the sample genome: decoys and the EBV genome (chrEBV, NC_007605) added to human references
    function is_foreign(c) {
        return c ~ /_decoy$|EBV|^hs37d5$|^NC_007605(\.[0-9]+)?$/
    }
    # Unplaced, alternate, decoy, HLA and EBV sequences, left out of the per-contig table. GenBank accessions
    # (GL000192.1, KI270706.1, JH584304.1) are unplaced contigs only in references that also have chromosomes:
    # a scaffold-level assembly can name all its scaffolds that way.
    function is_minor(c) {
        return is_foreign(c) || c ~ /(_alt|_random|_fix|_hap[0-9]*)$|^chrUn|^HLA-/ ||
            (has_chrom_names && c ~ /^(GL|KI|JH|KN|KZ|KQ|ML|MU)[0-9]+\.[0-9]+$/)
    }
    # Chromosome-like names: numbered, sex, arm-named (Drosophila 2L) or Roman-numbered (yeast, worm)
    function is_chrom_name(c) {
        return tolower(c) ~ /^(chr)?([0-9]+|[0-9]+[lr]|[ivx]+|[xywz])$/
    }
    # In the genome-wide figures: chromosomes and their unplaced or alternate contigs, or every contig if the
    # reference has no chromosome-like names, less organelles, decoys, EBV and the -c contigs
    function is_genome(c) {
        if (is_organelle(c) || is_foreign(c) || (c in control)) return 0
        return has_chrom_names ? (is_chrom_name(c) || is_minor(c)) : 1
    }
    function roman(s,    i, v, next_v, total) {
        total = 0
        for (i = 1; i <= length(s); i++) {
            v = roman_value[substr(s, i, 1)]; next_v = roman_value[substr(s, i + 1, 1)]
            total += (v < next_v) ? -v : v
        }
        return total
    }
    # Natural order: numbered and Roman-numbered chromosomes (2L before 2R), X, Y, organelles, then the rest by name
    function sort_key(c,    b) {
        b = c; sub(/^[Cc]hr/, "", b)
        if (b ~ /^[0-9]+$/) return sprintf("0%06d", b + 0)
        if (b ~ /^[0-9]+[LR]$/) return sprintf("0%06d%s", substr(b, 1, length(b) - 1) + 0, substr(b, length(b)))
        if (b == "X") return "1X"
        if (b == "Y") return "2Y"
        if (b ~ /^[IVX]+$/) return sprintf("0%06d", roman(b))
        if (is_organelle(c)) return "3" c
        return "4" c
    }
    # Order for the top contigs: most covered CpGs first, then by name so that ties are deterministic
    function before(a, b) {
        return contig_sites[a] > contig_sites[b] || (contig_sites[a] == contig_sites[b] && a < b)
    }
    BEGIN {
        sample = ENVIRON["sample"]; name = ENVIRON["name"]; label = ENVIRON["label"]
        n_controls = split(ENVIRON["control_contigs"], control_list, ",")
        for (i = 1; i <= n_controls; i++) if (control_list[i] != "") control[control_list[i]] = 1
        roman_value["I"] = 1; roman_value["V"] = 5; roman_value["X"] = 10; roman_value[""] = 0
        max_contigs = 50
        max_logged = 20
        for (i = 0; i < 10; i++) hist[i] = 0
        if (format == "bedgraph") { c_mod = 5; c_unmod = 6; source = "input"; ref = "" }
        else { source = "BED"; ref = "reference " }
    }
    format == "bedgraph" && /^#context\t/ { context = $2; next }
    format == "bedgraph" && /^track / { next }
    /^#/ {
        if (!rows && !c_cpg) {
            sub(/^#/, "")
            for (i = 1; i <= NF; i++) col[$i] = i
            c_unmod = col["unmod"]; c_mod = col["mod"]; c_cpg = col["cpg"]
        }
        next
    }
    format == "rastair2" {
        rows++
        if (!(c_unmod && c_mod && c_cpg)) {
            print "methylation_summary.sh: no unmod, mod and cpg columns in the BED header of " sample > "/dev/stderr"
            failed = 1; exit 1
        }
        if ($c_cpg != "REF") next
    }
    # CHG and CHH: only pooled per contig, for the genome-wide General Statistics columns
    format == "bedgraph" && context != "CpG" {
        seen[$1] = 1
        ctx_mod[context, $1] += $c_mod; ctx_depth[context, $1] += $c_mod + $c_unmod
        next
    }
    {
        seen[$1] = 1
        unmod = $c_unmod; mod = $c_mod; depth = unmod + mod
        if (depth == 0) next
        sites++
        contig_mod[$1] += mod; contig_depth[$1] += depth; contig_sites[$1]++
        if (depth >= min_depth) {
            bin = int(10 * mod / depth); if (bin > 9) bin = 9
            contig_hist[$1, bin]++
        }
    }
    END {
        if (failed) exit 1
        # No covered reference CpG: write nothing rather than a 0% methylation row in General Statistics
        if (sites == 0) {
            print "methylation_summary.sh: no covered " ref "CpGs for " sample "; no MultiQC rows written" > "/dev/stderr"
            exit 0
        }
        for (c in control) {
            if (!(c in seen)) print "methylation_summary.sh: " sample ": control contig " c " is not in the " source > "/dev/stderr"
        }
        # Decided from every CpG in the input, covered or not, so coverage does not change the choice
        has_chrom_names = 0
        for (c in seen) {
            if (!is_organelle(c) && !is_foreign(c) && !(c in control) && is_chrom_name(c)) { has_chrom_names = 1; break }
        }

        # Genome-wide figures, from the per-contig totals; log the contigs used, the first max_logged in natural order
        sites = 0; mod_sum = 0; depth_sum = 0; n_genome = 0; n_logged = 0
        for (c in contig_depth) {
            if (!is_genome(c)) continue
            n_genome++
            sites += contig_sites[c]; mod_sum += contig_mod[c]; depth_sum += contig_depth[c]
            for (i = 0; i < 10; i++) hist[i] += contig_hist[c, i]
            k = sort_key(c)
            if (n_logged == max_logged && k >= logged_key[n_logged]) continue
            if (n_logged < max_logged) n_logged++
            for (j = n_logged - 1; j >= 1 && k < logged_key[j]; j--) { logged[j + 1] = logged[j]; logged_key[j + 1] = logged_key[j] }
            logged[j + 1] = c; logged_key[j + 1] = k
        }
        if (n_genome == 0) {
            print "methylation_summary.sh: " sample ": no covered CpGs outside organelles, decoys, EBV and control contigs; no genome-wide figures written" > "/dev/stderr"
        }
        else {
            names = logged[1]
            for (i = 2; i <= n_logged; i++) names = names ", " logged[i]
            print "methylation_summary.sh: " sample ": genome-wide figures from " n_genome " contig(s): " names (n_genome > n_logged ? ", ..." : "") > "/dev/stderr"
            for (k in ctx_depth) {
                split(k, parts, SUBSEP)
                if (is_genome(parts[2])) { other_mod[parts[1]] += ctx_mod[k]; other_depth[parts[1]] += ctx_depth[k] }
            }
            write_genome_wide()
        }
        write_contig_table()
    }
    # General Statistics columns and the methylation distribution, from the genome-wide totals
    function write_genome_wide(    gs, hf, header, line, i, x, contexts, ctx) {
        gs = sample "." name "_generalstats_mqc.tsv"
        print "# id: \"" name "_generalstats\"" > gs
        print "# plot_type: \"generalstats\"" > gs
        print "# headers:" > gs
        print "#     " name "_cpg_meth:" > gs
        print "#         title: \"CpG meth\"" > gs
        print "#         namespace: \"" label "\"" > gs
        print "#         description: \"CpG methylation pooled over reads, " ref "CpGs on the genome, without organelles and spike-ins (" label ")\"" > gs
        print "#         suffix: \"%\"" > gs
        print "#         min: 0" > gs
        print "#         max: 100" > gs
        print "#         scale: \"RdYlBu-rev\"" > gs
        print "#         format: \"{:,.1f}\"" > gs
        print "#     " name "_cpg_sites:" > gs
        print "#         title: \"CpGs\"" > gs
        print "#         namespace: \"" label "\"" > gs
        print "#         description: \"" (ref ? "Reference CpG positions (per strand)" : "CpG positions") " with at least one informative read\"" > gs
        print "#         scale: \"Greens\"" > gs
        print "#         format: \"{:,.0f}\"" > gs
        print "#     " name "_cpg_depth:" > gs
        print "#         title: \"CpG depth\"" > gs
        print "#         namespace: \"" label "\"" > gs
        print "#         description: \"Mean informative reads per covered CpG position\"" > gs
        print "#         suffix: \"X\"" > gs
        print "#         scale: \"Blues\"" > gs
        print "#         format: \"{:,.1f}\"" > gs
        header = "Sample\t" name "_cpg_meth\t" name "_cpg_sites\t" name "_cpg_depth"
        line = sprintf("%s\t%.2f\t%d\t%.2f", sample, (depth_sum ? 100 * mod_sum / depth_sum : 0), sites, (sites ? depth_sum / sites : 0))
        split("CHG CHH", contexts, " ")
        for (x = 1; x <= 2; x++) {
            ctx = contexts[x]
            if (!other_depth[ctx]) continue
            print "#     " name "_" tolower(ctx) "_meth:" > gs
            print "#         title: \"" ctx " meth\"" > gs
            print "#         namespace: \"" label "\"" > gs
            print "#         description: \"" ctx " methylation pooled over reads on the genome, without organelles and spike-ins (" label ")\"" > gs
            print "#         suffix: \"%\"" > gs
            print "#         min: 0" > gs
            print "#         max: 100" > gs
            print "#         scale: \"RdYlBu-rev\"" > gs
            print "#         format: \"{:,.1f}\"" > gs
            header = header "\t" name "_" tolower(ctx) "_meth"
            line = line sprintf("\t%.2f", 100 * other_mod[ctx] / other_depth[ctx])
        }
        print header > gs
        print line > gs

        hf = sample "." name "_beta_hist_mqc.tsv"
        print "# id: \"" name "_beta_hist\"" > hf
        print "# section_name: \"" label " CpG methylation distribution\"" > hf
        print "# description: \"" (ref ? "Reference CpG" : "CpG") " positions on the genome (without organelles and spike-ins) in each methylation bin, for positions with at least " min_depth " informative reads.\"" > hf
        print "# plot_type: \"bargraph\"" > hf
        print "# pconfig:" > hf
        print "#     id: \"" name "_beta_hist_plot\"" > hf
        print "#     title: \"" label ": CpG methylation distribution\"" > hf
        print "#     ylab: \"CpG positions\"" > hf
        print "#     cpswitch_c_active: false" > hf
        header = "Sample"; line = sample
        for (i = 0; i < 10; i++) {
            header = header "\t" (10 * i) "-" (10 * i + 10) "%"
            line = line "\t" hist[i]
        }
        print header > hf
        print line > hf
    }
    # Per-contig table; MultiQC rejects a table without columns, so skip it when no contig qualifies.
    # Columns in natural order (see sort_key). Draft assemblies can have thousands of contigs: keep the
    # max_contigs with the most covered CpGs, in a list sorted as by before(), so each contig costs at
    # most max_contigs comparisons.
    function write_contig_table(    n_contigs, c, i, j, k, header, line, cf) {
        n_contigs = 0
        for (c in contig_depth) {
            if (is_minor(c) || contig_sites[c] < 10) continue
            if (n_contigs == max_contigs && !before(c, keep[n_contigs])) continue
            if (n_contigs < max_contigs) n_contigs++
            for (j = n_contigs - 1; j >= 1 && before(c, keep[j]); j--) keep[j + 1] = keep[j]
            keep[j + 1] = c
        }
        for (i = 1; i <= n_contigs; i++) key[i] = sort_key(keep[i])
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
            cf = sample "." name "_contigs_mqc.tsv"
            print "# id: \"" name "_contigs\"" > cf
            print "# section_name: \"" label " CpG methylation per contig\"" > cf
            print "# description: \"CpG methylation (%) pooled over reads for each contig, including any methylated or unmethylated spike-in controls. Contigs with fewer than 10 covered CpG positions, and unplaced, alternate, decoy, HLA and EBV sequences, are left out; at most " max_contigs " contigs (those with the most covered CpG positions) are shown.\"" > cf
            print "# plot_type: \"table\"" > cf
            print "# pconfig:" > cf
            print "#     id: \"" name "_contigs_table\"" > cf
            print "#     title: \"" label ": CpG methylation per contig\"" > cf
            print "#     min: 0" > cf
            print "#     max: 100" > cf
            print "#     format: \"{:,.1f}\"" > cf
            print header > cf
            print line > cf
        }
    }
'
