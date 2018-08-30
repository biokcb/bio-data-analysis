#/usr/bin/env python

###############################################
#
# trimming and tailing analysis functions
# 
# must have samtools and tophat installed
#
#
###############################################
import sys
import os.path
from lib.cmd_utils import run_cmd
from lib.file_util import grabLine, addLine, headerRep

# write accepted hits table
def appendTrimSuccess(outFile, loopNum, tophatFolder):
    """
    Append the successfully trimmed nucleotides
    to the header of a mapped sequence if the 
    sequence was mapped.
    """
    run_cmd("samtools bam2fq ./" + tophatFolder + "/accepted_hits.bam > mappedTEMP.fq")
    with open("mappedTemp.fq") as f:
        size=sum(1 for _ in f)
    f.close()
    with open("mappedTemp.fq") as f:
        fastq = f.readlines()
        lineNum = (0,1)
        countSeq=0
        while lineNum[0] < size:
            seqHeader, seq = grabLine(fastq, lineNum)
            countSeq += 1
            if loopNum > 1:
                seqName, trimSeq = seqHeader.split(":")
            else: 
                seqName = seqHeader
                trimSeq = ''
            finalList = [seqName+'\t', seq+'\t', str(loopNum)+'\t', trimSeq+'\n']
            addLine(finalList, outFile)	
            lineNum = (lineNum[0] + 4,lineNum[1] + 4)
    run_cmd("rm mappedTEMP.fq")
    print("Processed " + str(countSeq) + " mapped sequences!")
    f.close()

def trimHits(inFile, outFile, loopNum, add=False):
    """ 
    Trim a nucleotide of the end of a sequence
    and store the trimmed nucleotide in the 
    header for later. 
    """
    if add:
        with open(inFile) as f:
            size=sum(1 for _ in f)
        with open(inFile) as f:
            fastq = f.readlines()
            lineNum = (0,1,2,3)
            while lineNum[0] < size:
                seqHeader, seq, extraName, seqQuality = grabLine(fastq, lineNum)
                if loopNum > 1:
                    seqName, trimSeq = seqHeader.split(":")
                else: 
                    seqName = seqHeader
                    trimSeq = ''
                seqNew = seq[0:len(seq)-1]
                seqQNew = seqQuality[0:len(seqQuality)-1]
                finalList = [seqName+":"+seq[-1:]+trimSeq+'\n',seqNew+'\n', extraName+'\n', seqQNew+'\n']
                if len(seqNew) > 29:
                    addLine(finalList, "seqTemp.fq")
                lineNum = (lineNum[0] + 4,lineNum[1] + 4,lineNum[2] + 4,lineNum[3] + 4)
        f.close()
        run_cmd("rm " + inFile)
        run_cmd("mv seqTemp.fq " + outFile)
    else:
        run_cmd("samtools bam2fq ./" + inFile + "/unmapped.bam > unmappedTEMP.fq")
        with open("unmappedTemp.fq") as f:
            size=sum(1 for _ in f)
        with open("unmappedTemp.fq") as f:
            fastq = f.readlines()
            lineNum = (0,1,2,3)
            while lineNum[0] < size:
                seqHeader, seq, extraName, seqQuality = grabLine(fastq, lineNum)
                if loopNum > 1:
                    seqName, trimSeq = seqHeader.split(":")
                else: 
                    seqName = seqHeader
                    trimSeq = ''
                seqNew = seq[0:len(seq)-1]
                seqQNew = seqQuality[0:len(seqQuality)-1]
                finalList = [seqName+":"+seq[-1:]+trimSeq+'\n',seqNew+'\n', extraName+'\n', seqQNew+'\n']
                if len(seqNew) > 29:
                    addLine(finalList, outFile)
                lineNum = (lineNum[0] + 4,lineNum[1] + 4,lineNum[2] + 4,lineNum[3] + 4)
        f.close()
        run_cmd("rm unmappedTEMP.fq")


    
def main(fastqFileIN, outFile, tophatFolder):
    """ 
    main routine:
    aligns sequences, take the unaligned sequences,
    trims a nucleotide, realigns until the sequence
    is too short

    takes 3 arguments
    fastqFileIN = input fastq file to analyze
    outFile = output file name
    tophatFolder = name of tophat folder to use
    
    """
    fastqFile= fastqFileIN.split(".")[0] + "_fixed.fq"
    headerRep(fastqFileIN, fastqFile)

    # starting with sequences that did not map to WS230 (tophat was run 1x)
    loopNum = 1
    run_cmd("tophat --library-type fr-firststrand -N 0 -i 30 --read-gap-length 0 --max-insertion-length 0 --max-deletion-length 0 --read-edit-dist 0 -p 6 -o " + tophatFolder + str(loopNum) + " ./mcherry_ref/mcherry " + fastqFile)
    appendTrimSuccess(outFile, loopNum, tophatFolder+str(loopNum))
    trimHits(tophatFolder+str(loopNum),"sequences_"+str(loopNum)+".fq", loopNum, add=False)
    fileNum=loopNum

    # loop start
    # Take unaligned sequences, trim a nt off the end, then
    # realign the sequence. Repeat until the size is too small
    # to accurately match. Store trimmed sequence in header
    # to analyze later

    for loopNum in range(1,46):
        print('Now processing alignment iteration ' + str(loopNum) + '\n')
        fastqFile = "sequences_"+str(fileNum)+".fq"
        loopNum += 1
        run_cmd("tophat --no-coverage-search --library-type fr-firststrand -N 0 -i 30 --max-insertion-length 0 --read-gap-length 0 --max-deletion-length 0 --read-edit-dist 0 -p 6 -o " + tophatFolder + str(loopNum) + " ./mcherry_ref/mcherry " + fastqFile)
        if os.path.isfile(tophatFolder+str(loopNum)+"/accepted_hits.bam"):
            fileNum = loopNum
            appendTrimSuccess(outFile, loopNum, tophatFolder+str(loopNum))
            trimHits(tophatFolder+str(loopNum),"sequences_"+str(loopNum)+".fq", loopNum, add=False)
        else: 
            trimHits("sequences_"+str(fileNum)+".fq","sequences_"+str(fileNum)+".fq", loopNum, add=True)

    # remove temporary tophat files to avoid clogging things up
    run_cmd("rm -rf sequences_*")
    run_cmd("rm -rf " + tophatFolder+"*")

if __name__ == 'main':    
    fastqFileIN = sys.argv[1]
    outFile = sys.argv[2]
    tophatFolder = sys.argv[3]
    main(fastqFileIN, outFile, tophatFolder) 

