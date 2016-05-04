#!/usr/bin/env python

import argparse
import sys
import pysam

import logging
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

from aimhii.__about__ import (
    __author__, __copyright__, __email__, __license__, __summary__, __title__,
    __uri__, __version__)

def main():
    parser = argparse.ArgumentParser(description="Filter BAM/SAM file for eliminate all reads that are not putative fusions")
    parser.add_argument("IN_SAM_FILE", help="BAM/SAM file containing mapped reads")
    parser.add_argument("OUT_SAM_FILE", help="BAM/SAM file containing mapped reads")
    args = parser.parse_args()

    filter_chimeric_reads(args.IN_SAM_FILE, args.OUT_SAM_FILE)

    
def filter_chimeric_reads(in_sam_filename,out_sam_filename):
    junction_list = []
    in_samfile = pysam.Samfile(in_sam_filename)
    out_samfile = pysam.AlignmentFile(out_sam_filename, "wb",template=in_samfile)
    logger.debug("aread.qname, samfile.getrname(aread.tid), aread.cigarstring, aread.aend, aread.alen, aread.pos, aread.qend, aread.qlen, aread.qstart, tags['SA']")
    for aread in in_samfile:
        tags = dict(aread.tags)
        if "SA" in tags:
            if (aread.flag & 0x900) == 0 : #read is primary
                out_samfile.write(aread)
    out_samfile.close()

if __name__ == "__main__":
    main()
