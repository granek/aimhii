#!/usr/bin/env python

import argparse
import os
import re
import subprocess
import tempfile, shutil, atexit
import sys
import signal
import fileinput
import extract_chimeras

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def main():
    parser = argparse.ArgumentParser(description="Extract all reads in SAM_FILE that map to REF_NAME, or have pair that maps to it (including fusion matches)")
    parser.add_argument("REF_GENOME", help="FASTA file containing reference genome that reads will be mapped to.")
    parser.add_argument("INSERT_SEQ", help="FASTA file containing the sequence of the insertion DNA fragment.")
    parser.add_argument("ADAPTER_FILE", help="FASTA file containing sequences of adapters to be removed from reads.")
    parser.add_argument("FASTQ_FILE", nargs="+", help="One file for single-end, two for paired end.  Files may be gzipped.")
    parser.add_argument("-t", "--tmpdir", help="Directory for outputing intermediate files (default: %(default)s).")
    parser.add_argument("--threads", type=int, default=1, help="Number of parallel processed to run (default: %(default)s).")
    parser.add_argument("--minseed", type=int, default=10, help="Minimum seed length for bwa meme (default: %(default)s).")

    parser.add_argument("--outfile", metavar="TABLE_FILE", default=sys.stdout,
                        type=argparse.FileType('w'),
                        help="Output detailed table of insertions to %(metavar)s")
    parser.add_argument("--minreads", type=int, metavar="NUM_READS", default=2, 
                        help="Clusters consisting of less than %(metavar)s reads are ignored (default: %(default)s)")
    parser.add_argument("--maxgap", type=int, metavar="GAP_LENGTH", default=5000, 
                        help="Clusters separated by no more than %(metavar)s) bases are considered to be cluster doublets (default: %(default)s)")

    
    # parser.add_argument("--insert", type=file, help="The sequence of the inserted DNA fragment.  Used to determine how to organize junctions.")
    # parser.add_argument("REF_NAME", help="Extract all reads w")
    args = parser.parse_args()

    if len(args.FASTQ_FILE)>2:
        msg = "Maximum number of FASTQ files is 2 (for read pairs).  You supplied {0}".format(len(args.FASTQ_FILE))
        raise argparse.ArgumentTypeError(msg)

    if args.tmpdir:
        tmpdir = args.tmpdir
        if not os.path.exists(tmpdir):
            os.makedirs(tmpdir)
    else:
        tmpdir = tempfile.mkdtemp(suffix="aimhii")
        atexit.register(shutil.rmtree, tmpdir)

    concatentated_genome = os.path.join(tmpdir,"ref_and_insert.fa")
    concatenate_files(concatentated_genome, [args.REF_GENOME, args.INSERT_SEQ])
    
    ref_index = run_bwaindex(concatentated_genome)
    trimmed_fastqs = run_fastqmcf(args.ADAPTER_FILE, args.FASTQ_FILE,tmpdir)

    if len(args.FASTQ_FILE)==2:
        joined_filenames = run_fastqjoin(trimmed_fastqs,tmpdir)
        pair_bamname = run_bwamem(ref_index,joined_filenames[:2], os.path.join(tmpdir,"pair"),args.threads,args.minseed)
        joined_bamname = run_bwamem(ref_index,[joined_filenames[2]], os.path.join(tmpdir,"joined"),args.threads,args.minseed)
        merged_bamname = run_samtools_merge(os.path.join(tmpdir,"merged.bam"),pair_bamname,joined_bamname)
        final_bamname = merged_bamname
    elif len(args.FASTQ_FILE)==1:
        final_bamname = run_bwamem(ref_index,trimmed_fastqs, os.path.join(tmpdir,"single"),args.threads,args.minseed)

    run_samtools_index(final_bamname)

    extract_chimeras.run_analysis(final_bamname, args.minreads, args.maxgap, args.outfile, args.INSERT_SEQ, insertonly=True)


def run_fastqmcf(adapter_file, fastq_file_list, outdir):
    # fastq-mcf info/illumina_adapter1.fasta synthetic/syn_R1_001.fastq.gz synthetic/syn_R2_001.fastq.gz -o final_fastqs/syn_R1.tmp.gz -o final_fastqs/syn_R2.tmp.gz
    outfile_list = []
    for fastq_name in fastq_file_list:
        path_changed = os.path.join(outdir,(os.path.basename(fastq_name)))
        name_match = re.match("^(.+?)(\.fastq)?(\.gz)?$",path_changed,re.IGNORECASE)
        logger.debug(name_match.groups())
        # outname = re.sub("(.*)(\.fastq)?(\.gz)?",r'\1.trim\2\3',path_changed, re.IGNORECASE)
        outname = "".join([name_match.group(1), ".trim"] + [group for group in name_match.groups()[1:] if group != None])
        outfile_list.append(outname)
    if len(fastq_file_list) == 1:
        cmd = ("fastq-mcf", adapter_file, fastq_file_list[0], "-o", outfile_list[0])
    elif len(fastq_file_list) == 2:
        cmd = ("fastq-mcf", adapter_file, fastq_file_list[0], fastq_file_list[1], "-o",outfile_list[0], "-o",outfile_list[1])
    else:
        msg = "Maximum number of FASTQ files is 2 (for read pairs)."
        raise argparse.ArgumentTypeError(msg)
    print >>sys.stderr, "-"*20 + "TRIMMING READS" + "-"*20
    logger.debug("-"*20 + "TRIMMING READS" + "-"*20)
    logger.debug(cmd)
    # See https://blog.nelhage.com/2010/02/a-very-subtle-bug/ for explanation of "preexec_fn" below
    subprocess.check_output(cmd, preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    return outfile_list

##--------------
## RUN fastq-join
##--------------
def run_fastqjoin(fastq_file_list, outdir):
    # fastq-join final_fastqs/syn_R1.trim.fastq.gz final_fastqs/syn_R2.trim.fastq.gz -o final_fastqs/syn.%.tmp.gz
    # mv final_fastqs/syn.un1.tmp.gz final_fastqs/syn.un1.fastq.gz
    # mv final_fastqs/syn.un2.tmp.gz final_fastqs/syn.un2.fastq.gz
    # mv final_fastqs/syn.join.tmp.gz final_fastqs/syn.join.fastq.gz
    outfile_template = os.path.join(outdir, "%.fastq.gz")
    cmd = ("fastq-join", fastq_file_list[0], fastq_file_list[1], "-o", outfile_template)
    print >>sys.stderr, "-"*20 + "JOINING READS" + "-"*20
    logger.debug("-"*20 + "JOINING READS" + "-"*20)
    logger.debug(cmd)
    # See https://blog.nelhage.com/2010/02/a-very-subtle-bug/ for explanation of "preexec_fn" below
    subprocess.check_output(cmd, preexec_fn=lambda: signal.signal(signal.SIGPIPE, signal.SIG_DFL))
    return [outfile_template.replace("%",middle) for middle in ("un1", "un2", "join")]

def run_bwaindex(reference_genome):
    # bwa index -p bwa_bug_test/synthetic/syn_concat bwa_bug_test/synthetic/syn_concat.fa 
    index_prefix = os.path.splitext(reference_genome)[0]
    cmd = ("bwa", "index", "-p", index_prefix, reference_genome)
    print >>sys.stderr, "-"*20 + "INDEXING GENOME" + "-"*20
    logger.debug("-"*20 + "INDEXING GENOME" + "-"*20)
    logger.debug(cmd)
    subprocess.check_output(cmd)
    return index_prefix

def run_bwamem(genome_index, fastq_list,outfile_root,numthreads=1,minseed=10):
    # bwa mem -t 12 -k 10 synthetic/syn_concat final_fastqs/syn.un1.fastq.gz final_fastqs/syn.un2.fastq.gz | samtools view -Sb - > bwa/syn.pair_TMP
    # bwa mem -t 12 -k 10 synthetic/syn_concat final_fastqs/syn.join.fastq.gz | samtools view -Sb - > bwa/syn.join_TMP

    bwa_mem_cmd = ["bwa", "mem", "-t", str(numthreads), "-k", str(minseed), genome_index] + fastq_list
    print >>sys.stderr, "-"*20 + "RUNNING BWA" + "-"*20
    logger.debug("-"*20 + "RUNNING BWA" + "-"*20)
    logger.debug(bwa_mem_cmd)
    p_bwa_mem = subprocess.Popen(bwa_mem_cmd, stdout=subprocess.PIPE) #,stderr=errhandle)

    sam_view_cmd = ["samtools", "view", "-Su", "-"]
    p_sam_view = subprocess.Popen(sam_view_cmd, stdin=p_bwa_mem.stdout, stdout=subprocess.PIPE) #,stderr=errhandle)

    sam_sort_cmd = ["samtools", "sort", "-",outfile_root]
    p_sam_sort = subprocess.Popen(sam_sort_cmd, stdin=p_sam_view.stdout) #,stderr=errhandle)

    p_bwa_mem.stdout.close()
    p_sam_view.stdout.close()
    p_sam_sort.communicate()[0]
    return outfile_root+".bam"

##--------------------------------
## MERGE pair and joined BAM Files
##--------------------------------
# mkdir -p results
def run_samtools_merge(mergebam, inbam1, inbam2):
    # samtools merge bwa/syn.merge.bam bwa/syn.pair.bam bwa/syn.join.bam

    merge_cmd = ["samtools", "merge", "-f", mergebam, inbam1, inbam2]
    print >>sys.stderr, "-"*20 + "MERGING BAMS" + "-"*20
    logger.debug("-"*20 + "MERGING BAMS" + "-"*20)
    logger.debug(merge_cmd)
    subprocess.check_output(merge_cmd)
    return mergebam

def run_samtools_index(inbam):
    # samtools merge bwa/syn.merge.bam bwa/syn.pair.bam bwa/syn.join.bam
    index_cmd = ["samtools", "index", inbam]
    print >>sys.stderr, "-"*20 + "INDEXING BAM" + "-"*20
    logger.debug("-"*20 + "INDEXING BAM" + "-"*20)
    logger.debug(index_cmd)
    subprocess.check_output(index_cmd)
    return inbam

##--------------------------------
## Concatenate Files
##--------------------------------
def concatenate_files(outfilename,filelist):
    with open(outfilename, 'w') as fout:
        for line in fileinput.input(filelist):
            fout.write(line)

# make -f bwamem_chimeras.mk all SYNTHETIC="TRUE"  --dry-run

##--------------
## RUN fastq-mcf
##--------------
# fastq-mcf info/illumina_adapter1.fasta synthetic/syn_R1_001.fastq.gz synthetic/syn_R2_001.fastq.gz -o final_fastqs/syn_R1.tmp.gz -o final_fastqs/syn_R2.tmp.gz
# mv final_fastqs/syn_R1.tmp.gz final_fastqs/syn_R1.trim.fastq.gz
# mv final_fastqs/syn_R2.tmp.gz final_fastqs/syn_R2.trim.fastq.gz
# mkdir -p final_fastqs

##--------------
## RUN fastq-join
##--------------
# fastq-join final_fastqs/syn_R1.trim.fastq.gz final_fastqs/syn_R2.trim.fastq.gz -o final_fastqs/syn.%.tmp.gz
# mv final_fastqs/syn.un1.tmp.gz final_fastqs/syn.un1.fastq.gz
# mv final_fastqs/syn.un2.tmp.gz final_fastqs/syn.un2.fastq.gz
# mv final_fastqs/syn.join.tmp.gz final_fastqs/syn.join.fastq.gz

##------------------------
## INDEX reference genome
##------------------------
# bwa index synthetic/syn_concat.fa -p synthetic/syn_concat

##--------------
## RUN bwa mem
##--------------
# mkdir -p bwa
# bwa mem -t 12 -k 10 synthetic/syn_concat final_fastqs/syn.join.fastq.gz | samtools view -Sb - > bwa/syn.join_TMP
# mv bwa/syn.join_TMP bwa/syn.join.unsrtbam

# bwa mem -t 12 -k 10 synthetic/syn_concat final_fastqs/syn.un1.fastq.gz final_fastqs/syn.un2.fastq.gz | samtools view -Sb - > bwa/syn.pair_TMP
# mv bwa/syn.pair_TMP bwa/syn.pair.unsrtbam
# mkdir -p bwa

##--------------
## SORT BAM Files
##--------------
# samtools sort bwa/syn.join.unsrtbam  bwa/syn.join_TMP
# mv bwa/syn.join_TMP.bam bwa/syn.join.bam
# mkdir -p bwa
# samtools sort bwa/syn.pair.unsrtbam  bwa/syn.pair_TMP
# mv bwa/syn.pair_TMP.bam bwa/syn.pair.bam

##--------------------------------
## MERGE pair and joined BAM Files
##--------------------------------
# samtools merge bwa/syn.merge.bam bwa/syn.pair.bam bwa/syn.join.bam
# mkdir -p results


##--------------------------------
## Run extract_chimeras.py
##--------------------------------
# python2.7 scripts/extract_chimeras.py bwa/syn.merge.bam --junction results/syn.merge.junc.tmp --fusionreads results/syn.merge.fusread.tmp --insert synthetic/syn_insert.fa
# mv results/syn.merge.junc.tmp results/syn.merge.junc
# mv results/syn.merge.fusread.tmp results/syn.merge.fusread




if __name__ == "__main__":
    main()
