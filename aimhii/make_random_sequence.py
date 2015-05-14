#!/usr/bin/env python

import sys
import argparse
import numpy
import string

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from aimhii.__about__ import (
    __author__, __copyright__, __email__, __license__, __summary__, __title__,
    __uri__, __version__)


BASES = ("A","C", "G", "T")
# UPPER_BASES = letters.upper().split("")
HEAD_BASE = "A"
TAIL_BASE = "G"


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("--prefix", metavar="PREFIX", help="Use %(metavar)s as begining of sequence name", default="seq")
    parser.add_argument("--pad", metavar="PADLEN",type=int,default=0,
                        help="Pad beginging of sequence with %(metavar)s '{0}' , and end of sequence with %(metavar)s '{1}'".format(HEAD_BASE,TAIL_BASE))
    parser.add_argument("-n", "--numseqs", type=int, help="number of synthetic sequence to generate", default=2)
    parser.add_argument("-s", "--seqlen", type=int, help="length of synthetic sequence", default=2000)
    parser.add_argument("-o", "--output",type=argparse.FileType('w'),default=sys.stdout,metavar="OUTFILE")
    parser.add_argument("--alphabet",default=','.join(BASES),help="Alphabet to use (default: %(default)s)")
    parser.add_argument("--randseed",type=int,metavar="SEED",help="Seed the random number generator with %(metavar)s")
    parser.add_argument("--lower", action='store_true', help="Generate lowercase sequence")
    ## parser.add_argument("--upper", action='store_false', help="Generate uppercase sequence")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s v{0}".format(__version__))

    args = parser.parse_args()

    if args.randseed:
        numpy.random.seed(args.randseed)
    
    rec_list = []
    alphabet = map(string.strip,args.alphabet.split(","))
    if args.lower:
        bases = map(string.lower, alphabet)
        head = HEAD_BASE.lower()
        tail = TAIL_BASE.lower()
    else:
        bases = map(string.upper, alphabet)
        head = HEAD_BASE.upper()
        tail = TAIL_BASE.upper()
        
    for seq_num in range(args.numseqs):
        seq_string = ''.join(numpy.random.choice(bases, size=(args.seqlen - 2*args.pad), replace=True, p=None))
        
        seq_string = ''.join((head*args.pad, seq_string, tail*args.pad))
        id = args.prefix + str(seq_num)

        gen_rec=SeqRecord(Seq(seq_string,IUPAC.unambiguous_dna),id,'','')
        rec_list.append(gen_rec)
    SeqIO.write(rec_list, args.output, "fasta")

                                 



        
if __name__ == "__main__":
    main()
