#!/usr/bin/env bash
#
# Summarise per-CpG methylation calls as MultiQC custom content.
#
# Usage: methylation_summary.sh [-f rastair2|rastair|bedgraph] [-n name] [-l label] [-H] [-d min_depth] [-c control_contigs] [--] <sample> <file>...
#
#   -f  input format (default rastair2):
#         rastair2  one rastair version 2 CpG BED (`rastair call --bed`); columns found by name in the
#                   header (unmod, mod, cpg), and only reference CpGs (cpg == REF) counted
#         rastair   one rastair 0.8 call file (`*.rastair_call.txt`): the same columns without cpg, as it
#                   only reports reference CpGs
#         bedgraph  MethylDackel bedGraphs or Bismark coverage files (chrom, start, end, %, methylated,
#                   unmethylated), plain or gzipped. A file named *_CHG.bedGraph or *_CHH.bedGraph (.gz) holds that context;
#                   any other holds CpGs. Two kinds of file from a separate spike-in control alignment can be
#                   added: Bismark CX reports named *.controls.CX_report.txt[.gz], whose contigs are all
#                   controls, and Bismark alignment reports named controls.*_report.txt.
#   -n  tool name for the MultiQC IDs and column names: rastair2 (default), rastair, methyldackel or bismark
#   -l  label for one of several summaries per sample (e.g. a target panel): added to the file names, MultiQC
#       IDs and column names (lower case, with . and - as _) and to the section and column titles
#   -H  hide the General Statistics columns by default (they can still be shown in the report)
#   -d  minimum informative reads for a CpG to count in the methylation distribution (default 5)
#   -c  comma-separated spike-in or control contigs (e.g. lambda,pUC19) to leave out of the genome-wide figures
#
# Writes, for one sample:
# (with -l, <name> below is <name>_<label>)
#   <sample>.<name>_generalstats_mqc.tsv  CpG methylation, CpGs covered and mean depth over the genome, plus
#                                         CHG and CHH methylation if given (General Statistics)
#   <sample>.<name>_beta_hist_mqc.tsv     distribution of per-CpG methylation for CpGs with >= min_depth reads
#   <sample>.<name>_contigs_mqc.tsv       CpG methylation per contig (chromosomes and any spike-ins or controls),
#                                         plus CHG and CHH methylation for controls that have them
#   <sample>.<name>_controls_generalstats_mqc.tsv  reads on the controls, if given an alignment report
#
# Methylation is pooled over reads: sum(mod) / sum(mod + unmod). The genome-wide figures (General
# Statistics and the distribution) use chromosome-like contigs (chr1, 2, chrX, 2L, chrIV,
# chr1_KI270706v1_random, ...) if the input has any, otherwise every contig. Either way they leave out
# organelles (chrM, MT, chloroplast), decoys, EBV and the -c contigs, so that spike-ins and controls do
# not skew them. The contigs used are logged to stderr; if none of them has covered CpGs, the
# genome-wide files are not written. The per-contig table lists every contig with at least 10 covered
# CpGs (the 50 with the most, if there are more), except unplaced, alternate, decoy, HLA and EBV sequences.
# The -c contigs are always listed, first in line for the 50, if they have any covered CpG.

set -euo pipefail
# Locale-independent number formatting (some awks write 60,37 under a comma-decimal LC_NUMERIC)
export LC_ALL=C

usage() {
    echo "Usage: methylation_summary.sh [-f rastair2|rastair|bedgraph] [-n name] [-l label] [-H] [-d min_depth] [-c control_contigs] [--] <sample> <file>..." >&2
    exit 2
}

format=rastair2
name=rastair2
panel=""
export hidden=""
min_depth=5
# Passed through the environment: awk -v would interpret backslash escapes in the names
export control_contigs=""
while getopts "f:n:l:Hd:c:" opt; do
    case "$opt" in
        f) format=$OPTARG ;;
        n) name=$OPTARG ;;
        l) panel=$OPTARG ;;
        H) hidden=1 ;;
        d) min_depth=$OPTARG ;;
        c) control_contigs=$OPTARG ;;
        *) usage ;;
    esac
done
shift $((OPTIND - 1))
[ $# -ge 2 ] || usage
[[ $format == rastair2 || $format == rastair || $format == bedgraph ]] || { echo "methylation_summary.sh: -f must be rastair2, rastair or bedgraph, got '$format'" >&2; exit 2; }
[[ $format == bedgraph || $# -eq 2 ]] || { echo "methylation_summary.sh: -f $format takes one file" >&2; exit 2; }
[[ $name =~ ^[a-z0-9_]+$ ]] || { echo "methylation_summary.sh: -n must be lower-case letters, digits and underscores, got '$name'" >&2; exit 2; }
[[ -z $panel || $panel =~ ^[A-Za-z0-9._-]+$ ]] || { echo "methylation_summary.sh: -l must be letters, digits, '.', '_' and '-', got '$panel'" >&2; exit 2; }
[[ $min_depth =~ ^[0-9]+$ ]] || { echo "methylation_summary.sh: -d must be a whole number, got '$min_depth'" >&2; exit 2; }
case "$name" in
    methyldackel) label=MethylDackel ;;
    bismark) label=Bismark ;;
    *) label=$name ;;
esac
export title_suffix="" desc_label=$label
if [ -n "$panel" ]; then
    desc_label="$label, $panel"
    name="${name}_$(printf '%s' "$panel" | tr 'A-Z.-' 'a-z__')"
    label="$label ($panel)"
    title_suffix=" ($panel)"
fi
export sample=$1 name label
shift

# The bedgraph format has no context column: mark where each file starts, its kind and context.
# The leading newline ends a last line that has none.
for file in "$@"; do
    if [ "$format" = bedgraph ]; then
        case "$file" in
            *.controls.CX_report.txt | *.controls.CX_report.txt.gz)
                printf '\n#format\tcx\t%s\n' "$file" ;;
            controls.*_report.txt | */controls.*_report.txt)
                first=$(head -n 1 -- "$file")
                [[ $first == "Bismark report for: "* ]] || { echo "methylation_summary.sh: $file is not a Bismark alignment report" >&2; exit 2; }
                printf '\n#format\tbismark_report\t%s\n' "$file" ;;
            *)
                printf '\n#format\tcov\t%s\n' "$file"
                case "$file" in
                    *_CHG.bedGraph | *_CHG.bedGraph.gz) printf '#context\tCHG\n' ;;
                    *_CHH.bedGraph | *_CHH.bedGraph.gz) printf '#context\tCHH\n' ;;
                    *) printf '#context\tCpG\n' ;;
                esac ;;
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
    # Order for the top contigs: controls, then most covered CpGs, then by name so that ties are deterministic
    function before(a, b) {
        if ((a in control) != (b in control)) return (a in control)
        return contig_sites[a] > contig_sites[b] || (contig_sites[a] == contig_sites[b] && a < b)
    }
    function add_cpg(c, mod, unmod,    depth, bin) {
        depth = unmod + mod
        sites++
        contig_mod[c] += mod; contig_depth[c] += depth; contig_sites[c]++
        if (depth >= min_depth) {
            bin = int(10 * mod / depth); if (bin > 9) bin = 9
            contig_hist[c, bin]++
        }
    }
    BEGIN {
        sample = ENVIRON["sample"]; name = ENVIRON["name"]; label = ENVIRON["label"]
        title_suffix = ENVIRON["title_suffix"]; hidden = ENVIRON["hidden"]; desc_label = ENVIRON["desc_label"]
        n_controls = split(ENVIRON["control_contigs"], control_list, ",")
        for (i = 1; i <= n_controls; i++) if (control_list[i] != "") control[control_list[i]] = 1
        roman_value["I"] = 1; roman_value["V"] = 5; roman_value["X"] = 10; roman_value[""] = 0
        max_contigs = 50
        max_logged = 20
        for (i = 0; i < 10; i++) hist[i] = 0
        if (format == "bedgraph") { c_mod = 5; c_unmod = 6; source = "input"; ref = "" }
        else { source = (format == "rastair" ? "calls" : "BED"); ref = "reference " }
    }
    !NF { next }
    format == "bedgraph" && /^#format\t/ { fmt = $2; file = $3; fmt_seen[fmt] = 1; next }
    format == "bedgraph" && /^#context\t/ { context = $2; next }
    format == "bedgraph" && /^track / { next }
    # Bismark alignment report of the control alignment: reads analysed and aligned uniquely
    fmt == "bismark_report" {
        if ($1 == "Sequence pairs analysed in total:" || $1 == "Sequences analysed in total:") analysed += $2
        else if ($1 == "Number of paired-end alignments with a unique best hit:" || $1 == "Number of alignments with a unique best hit from the different alignments:") unique += $2
        next
    }
    # Bismark CX report (chrom, position, strand, methylated, unmethylated, context, trinucleotide): all controls
    fmt == "cx" {
        if ($4 !~ /^[0-9]+$/ || $5 !~ /^[0-9]+$/) {
            print "methylation_summary.sh: " file ": not a count in: " $0 > "/dev/stderr"
            failed = 1; exit 1
        }
        ctx = ($6 == "CG") ? "CpG" : $6
        if (ctx != "CpG" && ctx != "CHG" && ctx != "CHH") {
            print "methylation_summary.sh: " file ": unknown context " $6 > "/dev/stderr"
            failed = 1; exit 1
        }
        seen[$1] = 1; control[$1] = 1; control_input[$1] = 1; ctx_rows[ctx, $1] = 1
        if ($4 + $5 == 0) next
        if (ctx == "CpG") add_cpg($1, $4, $5)
        else { ctx_mod[ctx, $1] += $4; ctx_depth[ctx, $1] += $4 + $5 }
        next
    }
    format == "bedgraph" && !/^#/ {
        if ($c_mod !~ /^[0-9]+$/ || $c_unmod !~ /^[0-9]+$/) {
            print "methylation_summary.sh: " file ": not a count in: " $0 > "/dev/stderr"
            failed = 1; exit 1
        }
        main_input[$1] = 1
    }
    /^#/ {
        if (!rows && !c_unmod) {
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
    format == "rastair" {
        rows++
        if (!(c_unmod && c_mod)) {
            print "methylation_summary.sh: no unmod and mod columns in the header of " sample > "/dev/stderr"
            failed = 1; exit 1
        }
    }
    # CHG and CHH: only pooled per contig, for the genome-wide General Statistics columns
    format == "bedgraph" && context != "CpG" {
        seen[$1] = 1; ctx_rows[context, $1] = 1
        ctx_mod[context, $1] += $c_mod; ctx_depth[context, $1] += $c_mod + $c_unmod
        next
    }
    {
        seen[$1] = 1
        if ($c_unmod + $c_mod == 0) next
        add_cpg($1, $c_mod, $c_unmod)
    }
    END {
        if (failed) exit 1
        for (c in control_input) {
            if (c in main_input) {
                print "methylation_summary.sh: " sample ": control contig " c " also has calls from the main alignment; give the control a name that is not in the reference" > "/dev/stderr"
                exit 1
            }
        }
        if (analysed) write_controls_stats()
        else if (fmt_seen["bismark_report"]) print "methylation_summary.sh: " sample ": the control alignment report counts no reads; no % on controls written" > "/dev/stderr"
        # No covered reference CpG: write nothing rather than a 0% methylation row in General Statistics
        if (sites == 0) {
            print "methylation_summary.sh: no covered " ref "CpGs for " sample "; no methylation rows written" > "/dev/stderr"
            # Controls from a control alignment are still listed, with empty cells
            for (c in control_input) { write_contig_table(); break }
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
        print "#         title: \"CpG meth" title_suffix "\"" > gs
        print "#         namespace: \"" label "\"" > gs
        print "#         description: \"CpG methylation pooled over reads, " ref "CpGs on the genome, without organelles and spike-ins (" desc_label ")\"" > gs
        print "#         suffix: \"%\"" > gs
        print "#         min: 0" > gs
        print "#         max: 100" > gs
        print "#         scale: \"RdYlBu-rev\"" > gs
        print "#         format: \"{:,.1f}\"" > gs
        if (hidden) print "#         hidden: true" > gs
        print "#     " name "_cpg_sites:" > gs
        print "#         title: \"CpGs" title_suffix "\"" > gs
        print "#         namespace: \"" label "\"" > gs
        print "#         description: \"" (ref ? "Reference CpG positions (per strand)" : "CpG positions") " with at least one informative read\"" > gs
        print "#         scale: \"Greens\"" > gs
        print "#         format: \"{:,.0f}\"" > gs
        if (hidden) print "#         hidden: true" > gs
        print "#     " name "_cpg_depth:" > gs
        print "#         title: \"CpG depth" title_suffix "\"" > gs
        print "#         namespace: \"" label "\"" > gs
        print "#         description: \"Mean informative reads per covered CpG position\"" > gs
        print "#         suffix: \"X\"" > gs
        print "#         scale: \"Blues\"" > gs
        print "#         format: \"{:,.1f}\"" > gs
        if (hidden) print "#         hidden: true" > gs
        header = "Sample\t" name "_cpg_meth\t" name "_cpg_sites\t" name "_cpg_depth"
        line = sprintf("%s\t%.2f\t%d\t%.2f", sample, (depth_sum ? 100 * mod_sum / depth_sum : 0), sites, (sites ? depth_sum / sites : 0))
        split("CHG CHH", contexts, " ")
        for (x = 1; x <= 2; x++) {
            ctx = contexts[x]
            if (!other_depth[ctx]) continue
            print "#     " name "_" tolower(ctx) "_meth:" > gs
            print "#         title: \"" ctx " meth" title_suffix "\"" > gs
            print "#         namespace: \"" label "\"" > gs
            print "#         description: \"" ctx " methylation pooled over reads on the genome, without organelles and spike-ins (" desc_label ")\"" > gs
            print "#         suffix: \"%\"" > gs
            print "#         min: 0" > gs
            print "#         max: 100" > gs
            print "#         scale: \"RdYlBu-rev\"" > gs
            print "#         format: \"{:,.1f}\"" > gs
            if (hidden) print "#         hidden: true" > gs
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
    # General Statistics: of the first --control_reads reads, % with a unique alignment to the control FASTA
    function write_controls_stats(    cs) {
        cs = sample "." name "_controls_generalstats_mqc.tsv"
        print "# id: \"" name "_controls\"" > cs
        print "# plot_type: \"generalstats\"" > cs
        print "# headers:" > cs
        print "#     " name "_control_reads:" > cs
        print "#         title: \"% on controls" title_suffix "\"" > cs
        print "#         namespace: \"" label "\"" > cs
        print "#         description: \"Of the first --control_reads reads (read pairs for paired-end data), % with a unique alignment to the separate spike-in control FASTA (" desc_label ")\"" > cs
        print "#         suffix: \"%\"" > cs
        print "#         min: 0" > cs
        print "#         max: 100" > cs
        print "#         scale: \"Purples\"" > cs
        print "#         format: \"{:,.2f}\"" > cs
        if (hidden) print "#         hidden: true" > cs
        print "Sample\t" name "_control_reads" > cs
        printf "%s\t%.2f\n", sample, 100 * unique / analysed > cs
    }
    # Per-contig table; MultiQC rejects a table without columns, so skip it when no contig qualifies.
    # Columns in natural order (see sort_key). Draft assemblies can have thousands of contigs: keep the
    # max_contigs with the most covered CpGs, in a list sorted as by before(), so each contig costs at
    # most max_contigs comparisons.
    # Controls are always listed: an empty cell for a control without covered positions, and CHG and CHH columns
    # for a control with calls in those contexts.
    function write_contig_table(    n_contigs, c, i, j, k, x, ctx, contexts, header, line, cf, has_ctx) {
        n_contigs = 0
        for (c in seen) {
            if (!(c in control) && (is_minor(c) || !(c in contig_depth) || contig_sites[c] < 10)) continue
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
        split("CHG CHH", contexts, " ")
        for (i = 1; i <= n_contigs; i++) {
            c = keep[i]
            header = header "\t" c
            line = line "\t" ((c in contig_depth) ? sprintf("%.2f", 100 * contig_mod[c] / contig_depth[c]) : "")
            if (!(c in control)) continue
            for (x = 1; x <= 2; x++) {
                ctx = contexts[x]
                if (!((ctx, c) in ctx_rows)) continue
                has_ctx = 1
                header = header "\t" c " " ctx
                line = line "\t" (ctx_depth[ctx, c] ? sprintf("%.2f", 100 * ctx_mod[ctx, c] / ctx_depth[ctx, c]) : "")
            }
        }
        if (n_contigs > 0) {
            cf = sample "." name "_contigs_mqc.tsv"
            print "# id: \"" name "_contigs\"" > cf
            print "# section_name: \"" label " CpG methylation per contig\"" > cf
            print "# description: \"CpG methylation (%) pooled over reads for each contig, including any methylated or unmethylated spike-in controls, which are always shown. Other contigs with fewer than 10 covered CpG positions, and unplaced, alternate, decoy, HLA and EBV sequences, are left out; at most " max_contigs " contigs (the controls, then those with the most covered CpG positions) are shown." (has_ctx ? " For spike-in and control contigs, CHG and CHH methylation are shown too." : "") "\"" > cf
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
