#!/usr/bin/env python

from getopt import getopt
import sys
import os
from glob import glob
import gzip
from types import StringTypes


HELP_STRING = """Given a set of paired end read files, performs the trimming for RRBS with
diversity adapters.  Make sure the pattern names are in quotations.

    -1   pattern for forward read files (for example, "R1_1_edit3_6bp_BC??.fq")
    -2   pattern for reverse read files (for example, "R2_1_edit3_6bp_BC??.fq")

If a -2 parameter is not provided then it will be assumed that the data is single ended.
The trimmed files will be named for the original files + "_trimmed.fq".



VERSION: $Revision: 1.11 $

"""

# DEFAULTS
BASE_POSITIONS = ["A", "C", "G", "T", "N"]
NUM_BASES = len(BASE_POSITIONS)
SAVE_UNDIGESTED = False
SHOW_STATS = False  # When this is set to False, no report of the diversity will be printed.  This needs to be toggled in the code if you want to turn it on.
USED_DGG = False    # When this is set to False, the new version using RGG will be used for the analysis.  This version of RGG is more specific than the previously used DGG.  This needs to be toggled in the code if you want to turn it on.

# FUNCTION HELPERS
def getEmptyStats():
    d = {"filename":0,
            "nrBothRE":0,
            "nrFwdRE":0,
            "nrRevRE":0,
            "nrNoRE":0,
            "nrFwdREisC":0,
            "nrFwdREisT":0,
            "bdFwd4thBase":[0,0,0,0,0],
            "nrD0":0,
            "nrCGGCGG":0,
            "nrTGGTGG":0,
            "nrCGGTGG":0,
            "nrTGGCGG":0,
            "bdD1base1":[0,0,0,0,0],
            "bdD2base1":[0,0,0,0,0],
            "bdD2base2":[0,0,0,0,0],
            "bdD3base1":[0,0,0,0,0],
            "bdD3base2":[0,0,0,0,0],
            "bdD3base3":[0,0,0,0,0],
            "bdRev4thBase":[0,0,0,0,0],
            "nrRevD0":0,
            "bdRevD1base1":[0,0,0,0,0],
            "bdRevD2base1":[0,0,0,0,0],
            "bdRevD2base2":[0,0,0,0,0],
            "bdRevD3base1":[0,0,0,0,0],
            "bdRevD3base2":[0,0,0,0,0],
            "bdRevD3base3":[0,0,0,0,0]
            }

    return d

def FastqIterator(fh):
    """return an iterator of Records found in file handle, fh.
    """
    def readTotitle(fh, titleChar):
        """returns a tuple ([lines before the next title line], next tile line)
        """
        preLines = []
        while True:
            l = fh.readline().strip()
            if l.startswith(titleChar):
                return (preLines,l)
            elif l == '':
                return preLines,None
            else:
                preLines.append(l)

    if type(fh) in StringTypes:
        if (fh.endswith(".gz")):
            fh = gzip.open(fh,'r')
        else:
            fh = file(fh)

    preLines,nextTitleLine =readTotitle(fh,'@')

    while nextTitleLine != None:
        seqTitle = nextTitleLine[1:].rstrip()
        preLines,nextTitleLine=readTotitle(fh,'+')
        qualTitle = nextTitleLine[1:].rstrip()
        if len(qualTitle.strip()) > 0 and seqTitle != qualTitle:
            print seqTitle
            print preLines
            print qualTitle
            raise Exception("Error in parsing: @title sequence entry must be immediately followed by corresponding +title quality entry.")
        seqLines = preLines
        qualLines = []
        for i in range(len(seqLines)): # Quality characters should be the same length as the sequence
            qualLines.append( fh.readline().strip() )

        preLines,nextTitleLine=readTotitle(fh,'@')

        yield (seqTitle, ''.join(seqLines), ''.join(qualLines))

def trimOneRecord(fwdTitle, fwdSeq, fwdQual, revTitle, revSeq, revQual, outFwd, outRev, stats):


    if msp1 == True:
        fwdmPos = fwdSeq[:6].find("CGG")
        fwdnmPos = fwdSeq[:6].find("TGG")
    elif taq1 == True:
        fwdmPos = fwdSeq[:6].find("TGA")
        fwdnmPos = fwdSeq[:6].find("CGA")
    elif both_taq1_msp1 == True:
        fwdmPos = fwdSeq[:6].find("CGG")
        fwdnmPos = fwdSeq[:6].find("TGG")
        if fwdmPos == -1:
            fwdmPos = fwdSeq[:6].find("TGA")
        if fwdnmPos == -1:
            fwdnmPos = fwdSeq[:6].find("CGA")

    # figure out CGR position
    if revSeq == None:
        cgaPos = -1
    elif revSeq:
        cgaPos = revSeq[:6].find("CGA")
        if cgaPos == -1:
            if msp1 == True or both_taq1_msp1 == True:
                cgaPos = revSeq[:6].find("CGG")

    # if a read doesn't have a fwd and rev MSP site, end early without writing to output
    if fwdmPos == -1 and fwdnmPos == -1 and cgaPos == -1:
        stats["nrNoRE"] += 1
        if not SAVE_UNDIGESTED:
            return
    elif cgaPos == -1 and revSeq != None:
        stats["nrFwdRE"] += 1
        if not SAVE_UNDIGESTED:
            return
    elif cgaPos == -1 and revSeq == None:    ##except if it's a single read
        stats["nrFwdRE"] += 1
    elif fwdmPos == -1 and fwdnmPos == -1:
        stats["nrRevRE"] += 1
        if not SAVE_UNDIGESTED:
            return
    else:
        stats["nrBothRE"] += 1

    ## FORWARD READ
    # figure out ygg position
    s = fwdSeq[:6].upper()
    if s == "CGGCGG":
        yggPos = 0
        stats["nrCGGCGG"] += 1
    elif s == "TGGTGG":
        if USED_DGG:
            yggPos = 3
        else:
            yggPos = 0
        stats["nrTGGTGG"] += 1
    elif s == "CGGTGG":
        yggPos = 0
        stats["nrCGGTGG"] += 1
    elif s == "TGGCGG":
        if USED_DGG:
            yggPos = 3
        else:
            yggPos = 0
        stats["nrTGGCGG"] += 1
    else:
        yggPos = max(fwdmPos, fwdnmPos)

    if yggPos == 0:
        stats["nrD0"] += 1
    elif yggPos == 1:
        if fwdSeq[0] not in BASE_POSITIONS:
            print "Hmm... odd letter %s in %s at fwdSeq[0] when yggPos == 1" % (fwdSeq[1], fwdTitle)
        for i in range(NUM_BASES):
            if fwdSeq[0] == BASE_POSITIONS[i]:
                stats["bdD1base1"][i] += 1
    elif yggPos == 2:
        if fwdSeq[0] not in BASE_POSITIONS:
            print "Hmm... odd letter %s in %s at fwdSeq[0] when yggPos == 2" % (fwdSeq[1], fwdTitle)
        if fwdSeq[1] not in BASE_POSITIONS:
            print "Hmm... odd letter %s in %s at fwdSeq[1] when yggPos == 2" % (fwdSeq[1], fwdTitle)
        for i in range(NUM_BASES):
            if fwdSeq[0] == BASE_POSITIONS[i]:
                stats["bdD2base1"][i] += 1
        for i in range(NUM_BASES):
            if fwdSeq[1] == BASE_POSITIONS[i]:
                stats["bdD2base2"][i] += 1
    elif yggPos == 3:
        if fwdSeq[0] not in BASE_POSITIONS:
            print "Hmm... odd letter %s in %s at fwdSeq[0] when yggPos == 3" % (fwdSeq[0], fwdTitle)
        if fwdSeq[1] not in BASE_POSITIONS:
            print "Hmm... odd letter %s in %s at fwdSeq[1] when yggPos == 3" % (fwdSeq[1], fwdTitle)
        if fwdSeq[2] not in BASE_POSITIONS:
            print "Hmm... odd letter %s in %s at fwdSeq[2] when yggPos == 3" % (fwdSeq[2], fwdTitle)
        for i in range(NUM_BASES):
            if fwdSeq[0] == BASE_POSITIONS[i]:
                stats["bdD3base1"][i] += 1
        for i in range(NUM_BASES):
            if fwdSeq[1] == BASE_POSITIONS[i]:
                stats["bdD3base2"][i] += 1
        for i in range(NUM_BASES):
            if fwdSeq[2] == BASE_POSITIONS[i]:
                stats["bdD3base3"][i] += 1
    elif yggPos == -1:
        pass
    else:
        raise Exception("Impossible location for YGG: %s in %s" % (yggPos, fwdSeq))

    if yggPos >= 0:
        # track C vs. T in YGG
        if fwdSeq[yggPos] == "C":
            stats["nrFwdREisC"] += 1
        else:
            stats["nrFwdREisT"] += 1

        # track the stats on the 4th position
        for i in range(NUM_BASES):
            if fwdSeq[yggPos+3] == BASE_POSITIONS[i]:
                stats["bdFwd4thBase"][i] += 1

    if yggPos >= 0:  # if YGG position is found
        if revSeq == None:  # if this is a single read, trim only 5 bases from 3' end
            outFwd.write("@%s\n%s\n+\n%s\n" % (fwdTitle, fwdSeq[yggPos:-5], fwdQual[yggPos:-5]))  # trim off the diversity before YGG site and 5 bases from 3' end
        elif revSeq != None: # if this is a paired end read, trim 1 extra base from the 3' end, total of 6 bases so that alignment will be better for short fragments.
            outFwd.write("@%s\n%s\n+\n%s\n" % (fwdTitle, fwdSeq[yggPos:-6], fwdQual[yggPos:-6]))  # trim off the diversity before YGG site and 6 bases from the 3' end
    else: # if no YGG site is found
        if revSeq == None:  # if this is a single read, trim only 5 bases from 3' end
            outFwd.write("@%s\n%s\n+\n%s\n" % (fwdTitle, fwdSeq[:-5], fwdQual[:-5]))
        elif revSeq != None: # if this is a paired end read, trim 1 extra base from the 3' end, total of 6 bases so that alignment will be better for short fragments.
            outFwd.write("@%s\n%s\n+\n%s\n" % (fwdTitle, fwdSeq[:-6], fwdQual[:-6]))


    ##REVERSE READ
    if not revSeq:
        return

    # track the D0, D1, D2, D3 stats and the base composition of the trimmed bases
    if cgaPos == 0:
        stats["nrRevD0"] += 1
    elif cgaPos == 1:
        for i in range(NUM_BASES):
            if revSeq[0] == BASE_POSITIONS[i]:
                stats["bdRevD1base1"][i] += 1
    elif cgaPos == 2:
        for i in range(NUM_BASES):
            if revSeq[0] == BASE_POSITIONS[i]:
                stats["bdRevD2base1"][i] += 1
        for i in range(NUM_BASES):
            if revSeq[1] == BASE_POSITIONS[i]:
                stats["bdRevD2base2"][i] += 1
    elif cgaPos == 3:
        for i in range(NUM_BASES):
            if revSeq[0] == BASE_POSITIONS[i]:
                stats["bdRevD3base1"][i] += 1
        for i in range(NUM_BASES):
            if revSeq[1] == BASE_POSITIONS[i]:
                stats["bdRevD3base2"][i] += 1
        for i in range(NUM_BASES):
            if revSeq[2] == BASE_POSITIONS[i]:
                stats["bdRevD3base3"][i] += 1
    elif cgaPos == -1:
        # this will only happen if we keep going with reads not full digested
        pass
    else:
        raise Exception("Impossible location for CGA: %s in %s" % (cgaPos, revSeq))

    # track the stats on the 4th position
    if cgaPos >= 0:
        for i in range(NUM_BASES):
            if revSeq[cgaPos+3] == BASE_POSITIONS[i]:
                stats["bdRev4thBase"][i] += 1

    if cgaPos >= 0:  # if cga position is found
        outRev.write("@%s\n%s\n+\n%s\n" % (revTitle, revSeq[cgaPos+2:-6], revQual[cgaPos+2:-6])) # trim the cg from the 5' end and 6 bases from the 3' end
    else: # if no cga site is found
        outRev.write("@%s\n%s\n+\n%s\n" % (revTitle, revSeq[:-6], revQual[:-6])) # trim the cg from the 5' end and 6 bases from the 3' end




# retrieve the user parameters
fwdFilenames = None
revFilenames = None
statsFilename = None
usedDGG = USED_DGG
taq1 = False
both_taq1_msp1 = False
msp1 = True

try:
    optlist, args = getopt(sys.argv[1:], "h1:2:o:y:tb")
except:
    print "Error retrieving options"
    print ""
    print HELP_STRING
    sys.exit(1)

for (opt, opt_arg) in optlist:
    if opt == "-h":
        print ""
        print HELP_STRING
        sys.exit(1)
    elif opt == "-1":
        fwdFilenames = opt_arg
    elif opt == "-2":
        revFilenames = opt_arg
    elif opt == "-o":
        statsFilename = opt_arg
    elif opt == "-t":
        taq1 = True
        msp1 = False
    elif opt == "-b":
        both_taq1_msp1 = True
        msp1 = False

# check required parameters exist
if taq1 == True and both_taq1_msp1 == True:
    print "\nYou can only choose one restriction enzyme option.  If you want to select both enzymes, please enter the -b option only."
    print
    print HELP_STRING
    sys.exit(1)


if fwdFilenames == None:
    print "\nYou must provide both a fwd fastq filename or both fwd & rev filenames."
    print
    print HELP_STRING
    sys.exit(1)

if SHOW_STATS:
    if statsFilename == None:
        print "\nYou must provide an output filename"
        print
        print HELP_STRING
        sys.exit(1)

# do the actual work
fwdFiles = glob(fwdFilenames)
if revFilenames:
    revFiles = glob(revFilenames)
else:
    revFiles = []

if len(fwdFiles) == 0:
    print "There are no files that fit the fwd pattern '%s'" % (fwdFilenames)
    print
    print HELP_STRING
    sys.exit(1)

if len(revFiles) == 0:
    print "A file pattern wasn't provided for read 2.  Assuming the run is single ended."

print "Your files are:"
print "Fwd files:"
for f in fwdFiles:
    print "\t%s" % (f)

if revFilenames:
    print "Rev files:"
    for f in revFiles:
        print "\t%s" % (f)

if revFilenames:
    if len(fwdFiles) != len(revFiles):
        print "The fwd and rev files must have the same number of files."
        print
        print HELP_STRING
        sys.exit(1)

if SHOW_STATS:
    statsOut = open(statsFilename, "w")
    statsOut.write("\t".join(["Filename", "Both_RE", "Only_Fwd_RE", "Only_Rev_RE",
               "No_RE", "Fwd_YG(G/A)_is_C", "Fwd_YG(G/A)_is_T",
               "Fwd_4th_is_A", "Fwd_4th_is_C", "Fwd_4th_is_G", "Fwd_4th_is_T", "Fwd_4th_is_N",
               "Rev_4th_is_A", "Rev_4th_is_C", "Rev_4th_is_G", "Rev_4th_is_T", "Rev_4th_is_N",
               "Fwd_CGGCGG", "Fwd_TGGTGG", "Fwd_CGGTGG", "Fwd_TGGCGG",
               "D0_reads",
               "D1_base1_is_A", "D1_base1_is_C", "D1_base1_is_G", "D1_base1_is_T", "D1_base1_is_N",
               "D2_base1_is_A", "D2_base1_is_C", "D2_base1_is_G", "D2_base1_is_T", "D2_base1_is_N",
               "D2_base2_is_A", "D2_base2_is_C", "D2_base2_is_G", "D2_base2_is_T", "D2_base2_is_N",
               "D3_base1_is_A", "D3_base1_is_C", "D3_base1_is_G", "D3_base1_is_T", "D3_base1_is_N",
               "D3_base2_is_A", "D3_base2_is_C", "D3_base2_is_G", "D3_base2_is_T", "D3_base2_is_N",
               "D3_base3_is_A", "D3_base3_is_C", "D3_base3_is_G", "D3_base3_is_T", "D3_base3_is_N",
               "Rev_D0_reads",
               "Rev_D1_base1_is_A", "Rev_D1_base1_is_C", "Rev_D1_base1_is_G", "Rev_D1_base1_is_T", "Rev_D1_base1_is_N",
               "Rev_D2_base1_is_A", "Rev_D2_base1_is_C", "Rev_D2_base1_is_G", "Rev_D2_base1_is_T", "Rev_D2_base1_is_N",
               "Rev_D2_base2_is_A", "Rev_D2_base2_is_C", "Rev_D2_base2_is_G", "Rev_D2_base2_is_T", "Rev_D2_base2_is_N",
               "Rev_D3_base1_is_A", "Rev_D3_base1_is_C", "Rev_D3_base1_is_G", "Rev_D3_base1_is_T", "Rev_D3_base1_is_N",
               "Rev_D3_base2_is_A", "Rev_D3_base2_is_C", "Rev_D3_base2_is_G", "Rev_D3_base2_is_T", "Rev_D3_base2_is_N",
               "Rev_D3_base3_is_A", "Rev_D3_base3_is_C", "Rev_D3_base3_is_G", "Rev_D3_base3_is_T", "Rev_D3_base3_is_N"]))
    statsOut.write("\n")

for i in range(len(fwdFiles)):
    ##print "Working on pair fwd='%s' and rev='%s'" % (fwdFiles[i], revFiles[i])
    recCount = 1
    fwdIt = FastqIterator(fwdFiles[i])
    (fwdTitle, fwdSeq, fwdQual) = fwdIt.next()
    if revFilenames:
        revIt = FastqIterator(revFiles[i])
        (revTitle, revSeq, revQual) = revIt.next()
        revRoot = os.path.splitext(revFiles[i])[0]
        if (revFiles[i].endswith(".gz")):
            outRev = gzip.open(revRoot + "_trimmed.fq.gz", "wb")
        else:
            outRev = open(revRoot + "_trimmed.fq", "w")
    else:
        revIt = None
        revTitle = revSeq = revQual = None
        outRev = None

    if revFilenames and fwdTitle.split()[0] != revTitle.split()[0]:
        print "The fwd title and rev title don't match"
        print "fwd title: '%s'" % (fwdTitle.split()[0])
        print "rev title: '%s'" % (revTitle.split()[0])
        raise Exception("fwd and rev titles don't match.  Not paired end sequence.")

    fwdRoot = os.path.splitext(fwdFiles[i])[0]
    if (fwdFiles[i].endswith(".gz")):
        outFwd = gzip.open(fwdRoot + "_trimmed.fq.gz", "wb")
    else:
        outFwd = open(fwdRoot + "_trimmed.fq", "w")
    stats = getEmptyStats()
    stats["filename"] = fwdRoot
    statsKeys = stats.keys()
    statsKeys.sort()


    while True:

        # trim one record and add to output
        trimOneRecord(fwdTitle, fwdSeq, fwdQual, revTitle, revSeq, revQual, outFwd, outRev, stats)

        try:
            (fwdTitle, fwdSeq, fwdQual) = fwdIt.next()
            if revFilenames:
                (revTitle, revSeq, revQual) = revIt.next()
        except StopIteration:
            break

        recCount += 1

    statsList = [ stats["filename"], stats["nrBothRE"], stats["nrFwdRE"],
                  stats["nrRevRE"], stats["nrNoRE"], stats["nrFwdREisC"],
                  stats["nrFwdREisT"] ]
    statsList.extend(stats["bdFwd4thBase"])
    statsList.extend(stats["bdRev4thBase"])
    statsList.append(stats["nrCGGCGG"])
    statsList.append(stats["nrTGGTGG"])
    statsList.append(stats["nrCGGTGG"])
    statsList.append(stats["nrTGGCGG"])
    statsList.append(stats["nrD0"])
    statsList.extend(stats["bdD1base1"])
    statsList.extend(stats["bdD2base1"])
    statsList.extend(stats["bdD2base2"])
    statsList.extend(stats["bdD3base1"])
    statsList.extend(stats["bdD3base2"])
    statsList.extend(stats["bdD3base3"])
    statsList.append(stats["nrRevD0"])
    statsList.extend(stats["bdRevD1base1"])
    statsList.extend(stats["bdRevD2base1"])
    statsList.extend(stats["bdRevD2base2"])
    statsList.extend(stats["bdRevD3base1"])
    statsList.extend(stats["bdRevD3base2"])
    statsList.extend(stats["bdRevD3base3"])
    if SHOW_STATS:
        statsOut.write("\t".join([str(x) for x in statsList]))
        statsOut.write("\n")

    print "\tDone with this pair. there were %s records" % (recCount)
    d1 = sum(stats["bdD1base1"])
    d2 = sum(stats["bdD2base1"])
    d3 = sum(stats["bdD3base1"])
    if revFilenames:
        other = stats["nrNoRE"] + stats["nrRevRE"] + stats["nrFwdRE"]
    else:
        other = stats["nrNoRE"]
    print "\tFwd:  D0:%s  D1:%s  D2:%s  D3:%s other=%s (total=%s)" % (stats["nrD0"],
                    d1, d2, d3, other, stats["nrD0"]+d1+d2+d3+other)
