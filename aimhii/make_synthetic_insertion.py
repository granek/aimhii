#!/usr/bin/env python

import sys
import argparse
import random

## from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from aimhii.__about__ import (
    __author__, __copyright__, __email__, __license__, __summary__, __title__,
    __uri__, __version__)

def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("REFSEQ", type=file, help="The sequence to be insertionally mutated (FASTA format)")
    parser.add_argument("INSERTSEQ", type=file, help="The sequence to be inserted (FASTA format)")
    # parser.add_argument("-n", "--numreads", type=int, help="number of synthetic reads to generate", default=10)
    # parser.add_argument("-r", "--readlen", type=int, help="length of synthetic reads", default=250)
    parser.add_argument("--invert", action="store_true", help="Invert INSERTSEQ before inserting", default=False)
    parser.add_argument("--delete",type=int,help="Number of bases to delete at insertion site (default: %(default)s",default=0)
    parser.add_argument("-n", "--numchroms",type=int,help="Number of chromomsome to mutate",default=1)
    parser.add_argument("-p", "--position",type=int,help="Position in REFSEQ to insert INSERTSEQ",default=-1)
    parser.add_argument("-o", "--output",type=argparse.FileType('w'),default=sys.stdout,metavar="OUTFILE")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s v{0}".format(__version__))
    args = parser.parse_args()

    insert_rec = SeqIO.read(args.INSERTSEQ,"fasta")
    # ref_iter = SeqIO.parse(args.REFSEQ,"fasta")
    ref_dict = SeqIO.index(args.REFSEQ.name, "fasta")
    mutated_rec_i = random.sample(xrange(len(ref_dict)), args.numchroms)
    mutated_chrom_list = []
    print len(ref_dict)
    for i, cur_key in enumerate(ref_dict):
        if i in mutated_rec_i:
            cur_rec = ref_dict[cur_key]
            if args.invert:
                mutated_id = "_".join(map(str,(cur_rec.id,args.position,insert_rec.id,"inv", args.position+len(insert_rec),cur_rec.id)))
                insert_rec = insert_rec.reverse_complement()
            else:
                mutated_id = "_".join(map(str,(cur_rec.id,args.position,insert_rec.id,args.position+len(insert_rec),cur_rec.id)))
            mutated_chrom = generate_insertion(cur_rec,insert_rec,args.position,args.delete)
            mutated_chrom.description = ""
            mutated_chrom.id = mutated_id
            
            mutated_chrom_list.append(mutated_chrom)
        else:
            mutated_chrom_list.append(cur_rec)

    SeqIO.write(mutated_chrom_list, args.output, "fasta")

    # HERE HERE HERE HERE HERE HERE HERE HERE
    
    # gen_rec_list = []
    # quals = [40] * args.readlen
    # for i in range(args.numreads):
    #     a_frag_len = random.randint(0, args.readlen)
    #     if args.notchimera:
    #         a_frag_len = args.readlen
    #     a_frag, a_info = gen_random_frag(rec_list,a_frag_len)
    #     b_frag, b_info = gen_random_frag(rec_list,args.readlen-a_frag_len)
    #     gen_rec=SeqRecord(a_frag+b_frag,'{0[0]}{0[1]}_{0[2]}{0[3]}___{1[0]}{1[1]}_{1[2]}{1[3]}'.format(a_info,b_info),'','')
    #     gen_rec.letter_annotations["phred_quality"] =  quals
    #     gen_rec_list.append(gen_rec)
    # SeqIO.write(gen_rec_list, args.output, args.format)

    # quals = rec.letter_annotations["phred_quality"]
    # SeqIO.write(itertools.chain.from_iterable(seq_iter_list), output_handle, "fastq-sanger")


def generate_insertion(refseq,insert_seq,insert_pos,delete=0):
    if insert_pos < 0:
        insert_pos = random.randrange(len(refseq))
    mutated_seq = refseq[0:insert_pos] + insert_seq + refseq[insert_pos+delete:]
    # HERE HERE HERE HERE HERE HERE HERE HERE
        
    # source = random.choice(seq_list)
    # left = random.randint(0, len(source)- frag_len)
    # right = left+frag_len
    # strand = random.choice(("-","+"))
    # frag = source.seq[left:right]
    # info = source.id,left+1,right,strand
    # if strand == "-":
    #     frag = frag.reverse_complement()
        
    return mutated_seq

def gen_random_frag(seq_list,frag_len):
    source = random.choice(seq_list)
    left = random.randint(0, len(source)- frag_len)
    right = left+frag_len
    strand = random.choice(("-","+"))
    frag = source.seq[left:right]
    info = source.id,left+1,right,strand
    if strand == "-":
        frag = frag.reverse_complement()
        
    return frag,info
                                 



        
if __name__ == "__main__":
    main()
