#!/usr/bin/env bash
# Widen each target 1 bp at an edge that cuts a CpG; write widened_targets.bed and target_issues.txt (skipped/clipped)
# Usage: widen_targets_cpg.sh <targets.bed[.gz]> <genome.fa> <genome.fa.fai>

set -euo pipefail

if [ "$#" -ne 3 ]; then
    echo "Usage: $(basename "$0") <targets.bed[.gz]> <genome.fa> <genome.fa.fai>" >&2
    exit 2
fi
targets=$1
fasta=$2
fai=$3

# getfasta looks for <fasta>.fai and builds its own index when it is missing
if [ ! -e "${fasta}.fai" ]; then
    ln -s "$fai" "${fasta}.fai"
fi

# Gzipped by content, not by name
if gzip -cd "$targets" > /dev/null 2>&1; then
    read_targets() { gzip -cd "$targets"; }
else
    read_targets() { cat "$targets"; }
fi

# Usable targets as contig, start, end, id, contig length; issues as reason, contig, start, end
read_targets | awk 'BEGIN { FS = OFS = "\t"; issues = "target_issues.txt"; printf "" > issues }
    FNR == NR { len[$1] = $2; next }
    { sub(/\r$/, "") }
    /^[ \t]*$/ || /^#/ || /^(track|browser)([ \t]|$)/ { next }
    {
        if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/) { print "skipped: coordinates are not integers", $1, $2, $3 > issues; next }
        if (!($1 in len)) { print "skipped: contig not in the FASTA", $1, $2, $3 > issues; next }
        s = $2 + 0; e = $3 + 0; L = len[$1] + 0
        if (s >= e) { print "skipped: start is not before end", $1, $2, $3 > issues; next }
        if (s >= L) { print "skipped: start is past the contig end (" L ")", $1, $2, $3 > issues; next }
        if (e > L) { print "clipped: end is past the contig end (" L ")", $1, $2, $3 > issues; e = L }
        print $1, s, e, ++n, L
    }' "$fai" - > usable.bed

if [ ! -s usable.bed ]; then
    if [ -s target_issues.txt ]; then
        echo "ERROR: no target in '${targets}' fits the FASTA:" >&2
        head -n 20 target_issues.txt >&2
    else
        echo "ERROR: '${targets}' has no target lines" >&2
    fi
    exit 1
fi

# 2 bp windows across each edge: <id>s spans the start, <id>e the end (none at a contig's first base or end)
awk 'BEGIN { FS = OFS = "\t" }
    {
        if ($2 > 0) print $1, $2 - 1, $2 + 1, $4 "s"
        if ($3 < $5) print $1, $3 - 1, $3 + 1, $4 "e"
    }' usable.bed > windows.bed

if [ -s windows.bed ]; then
    bedtools getfasta -fi "$fasta" -bed windows.bed -bedOut > windows_seq.bed
else
    : > windows_seq.bed
fi

# getfasta skips a window it can't read with exit 0
n_windows=$(wc -l < windows.bed)
n_seq=$(wc -l < windows_seq.bed)
if [ "$n_windows" -ne "$n_seq" ]; then
    echo "ERROR: bedtools getfasta returned ${n_seq} of ${n_windows} edge windows" >&2
    exit 1
fi

# Window lookups are in target order, so read them in step with the targets
awk 'BEGIN { FS = OFS = "\t"; seqs = "windows_seq.bed" }
    function lookup(id,    line, f) {
        if ((getline line < seqs) <= 0) { print "ERROR: edge window " id " is missing" > "/dev/stderr"; exit 1 }
        split(line, f, "\t")
        if (f[4] != id) { print "ERROR: expected edge window " id ", found " f[4] > "/dev/stderr"; exit 1 }
        return toupper(f[5])
    }
    {
        s = $2; e = $3
        if ($2 > 0 && lookup($4 "s") == "CG") s = $2 - 1
        if ($3 < $5 && lookup($4 "e") == "CG") e = $3 + 1
        print $1, s, e
    }' usable.bed > widened_targets.bed
