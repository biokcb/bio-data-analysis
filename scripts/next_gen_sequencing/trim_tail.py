#!/usr/bin/env python
"""
A script for analyzing the trimming and tailing of mRNA sequencing data for signatures
of tailing or trimming in mRNA targeted by small RNA factors.

This script requires TopHat2 and samtools.
"""
import argparse
import os.path
import subprocess

def get_args():
    """
    Parse arguments from user input/commandline

    Requires that the user provide an input fastq file, an output file, and a tophat folder name
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input-fastq', metavar='FASTQ_FILE', required=True,
                        help='Input fastq file to process tails')
    parser.add_argument('-o', '--output-file', metavar='OUTPUT_FILE', required=True,
                        help='Output file name')
    parser.add_argument('-t', '--tophat-folder', metavar='TOPHAT_FOLDER', required=True,
                        help='Prefix for naming tophat folders')

    args = parser.parse_args()

    return args

def run_cmd(cmd):
    """
    Runs the requested command.

    Inputs:
        cmd: command to run

    Outputs:
        stdout: standard out from command
        stderr: standard error from command
    """
    print('Now running:')
    print(cmd)
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    output, error = process.communicate()

    return output, error

def read_fastq(fastq_file):
    """
    Reads the fastq file, assuming 4 lines per entry.

    Inputs:
        fastq_file: the open file handle

    Outputs:
        yields the four lines per sequence
    """
    while True:
        # Read in 4 lines at a time
        seq_header = f.readline().strip()
        seq = f.readline().strip()
        seq_comment = f.readline().strip()
        seq_quality = f.readline().strip()

        if not seq_quality:
            break

        yield seq_header, seq, seq_comment, seq_quality

def replace_header(fastq_file, output_file):
    """
    Replace header of a fastq file with a simpler name for compatibility with certain tools.

    Inputs:
        fastq_file = input fastq file to clean
        output_file = output file name to write

    Outputs:
        A fastq file named by 'output_file' with new headers:
        @seq1 .... @seqN
    """
    with open(fastq_file), open(output_file, 'a') as infile, outfile:
        seq_count = 0
        for seq_header, seq, seq_comment, seq_quality in read_fastq(infile):
            # Replace the header & write to the new file
            seq_header = '@seq_' + str(seq_count)
            final_list = '\n'.join([seq_header, seq, seq_comment, seq_quality])
            outfile.writelines(final_list)
            seq_count += 1

def append_trim_success(output_file, i, tophat_folder):
    """
    Append the successfully trimmed nucleotides to the header of a mapped sequence if the
    sequence was mapped.

    Inputs:
        output_file: The file to write appended hits to
        i: the loop number
        tophat_folder: the tophat folder to read aligned hits from
    """
    # Convert to fastq
    run_cmd("samtools bam2fq ./" + tophat_folder + "/accepted_hits.bam > mappedTEMP.fq")
    
    with open("mappedTemp.fq"), open(output_file, 'a') as mapped, outfile:
        seq_count = 0
        for seq_header, seq, seq_comment, seq_quality in read_fastq(infile):
            if i > 1:
                seq_name, trim_seq = seq_header.split(":")
            else:
                seq_name = seq_header
                trim_seq = ''
            final_list = [seq_name + '\t', seq + '\t', str(i) + '\t', trim_seq + '\n']
            outfile.writelines(final_list)
            seq_count += 1
    
    run_cmd("rm mappedTEMP.fq")
    print("Processed " + str(seq_count) + " mapped sequences!")

def trim_hits(in_file, output_file, i, add=False):
    """
    Trim a nucleotide of the end of a sequence and store the trimmed nucleotide in the
    header for later.

    Inputs:
        in_file: the input file to process
        output_file: the output file to write to
        i: the loop number
        add: boolean indicating whether or not it is being added to an existing fq
    """
    if add:
        with open(in_file), open("seqTemp.fq", 'a') as infile, outfile:
            for seq_header, seq, seq_comment, seq_quality in read_fastq(infile):
                if i > 1:
                    seq_name, trim_seq = seq_header.split(":")
                else:
                    seq_name = seq_header
                    trim_seq = ''
                new_seq = seq[0:len(seq)-1]
                new_seq_qual = seq_quality[0:len(seq_quality)-1]
                final_list = '\n'.join([seq_name + ":" + seq[-1:] + trim_seq,
                                        new_seq, seq_comment, new_seq_qual])
                if len(new_seq) > 29:
                    outfile.writelines(final_list)

        run_cmd("rm " + in_file)
        run_cmd("mv seqTemp.fq " + output_file)

    else:
        run_cmd("samtools bam2fq ./" + in_file + "/unmapped.bam > unmappedTEMP.fq")
        with open('unmappedTemp.fq'), open("seqTemp.fq", 'a') as infile, outfile:
            for seq_header, seq, seq_comment, seq_quality in read_fastq(infile):
                if i > 1:
                    seq_name, trim_seq = seq_header.split(":")
                else:
                    seq_name = seq_header
                    trim_seq = ''
                new_seq = seq[0:len(seq)-1]
                new_seq_qual = seq_quality[0:len(seq_quality)-1]
                final_list = '\n'.join([seq_name + ":" + seq[-1:] + trim_seq,
                                        new_seq, seq_comment, new_seq_qual])
                if len(new_seq) > 29:
                    outfile.writelines(final_list)
        run_cmd("rm unmappedTEMP.fq")

def main():
    """
    Main Routine:

    Aligns sequences, take the unaligned sequences, trims a nucleotide, realigns until the sequence
    is too short
    """
    args = get_args()

    # First replace the headers in the fastq file
    new_fastq_file = args.input_fastq.split(".")[0] + "_fixed.fq"
    replace_header(args.input_fastq, new_fastq_file)

    # Start with sequences that did not map to WS230 (tophat was run 1x)
    i = 1
    run_cmd("""
            tophat --library-type fr-firststrand 
            -N 0 
            -i 30 
            --read-gap-length 0 
            --max-insertion-length 0 
            --max-deletion-length 0 
            --read-edit-dist 0 
            -p 6 
            -o """ + args.tophat_folder + str(i)
            + " ./mcherry_ref/mcherry " + new_fastq_file)

    append_trim_success(args.output_file, i, args.tophat_folder + str(i))
    trim_hits(args.tophat_folder + str(i), "sequences_" + str(i) + ".fq", i, add=False)
    file_num = i

    # Take unaligned sequences, trim a nt off the end, then realign the sequence.
    # Repeat until the size is too small to accurately match. Store trimmed sequence
    # in header to analyze later
    for i in range(1, 30):
        print('Now processing alignment iteration ' + str(i) + '\n')
        new_fastq_file = "sequences_" + str(file_num) + ".fq"
        i += 1
        run_cmd("""
                tophat 
                --no-coverage-search 
                --library-type fr-firststrand 
                -N 0 
                -i 30 
                --max-insertion-length 0 
                --read-gap-length 0 
                --max-deletion-length 0 
                --read-edit-dist 0 
                -p 6 
                -o """ + args.tophat_folder + str(i)
                + " ./mcherry_ref/mcherry " + new_fastq_file)

        # Only append the trimmed hits if it aligned
        if os.path.isfile(args.tophat_folder + str(i) + "/accepted_hits.bam"):
            file_num = i
            append_trim_success(args.output_file, i, args.tophat_folder + str(i))
            trim_hits(args.tophat_folder+str(i), "sequences_" + str(i) + ".fq", i, add=False)
        else:
            trim_hits("sequences_" + str(file_num) + ".fq",
                      "sequences_" + str(file_num) + ".fq", i, add=True)

    # remove temporary tophat files to avoid clogging things up
    run_cmd("rm -rf sequences_*")
    run_cmd("rm -rf " + args.tophat_folder + "*")

if __name__ == '__main__':
    main()
