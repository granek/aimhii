import HTSeq
import argparse
import sys
import pysam
import re
import os

cigar_re = re.compile("(d+)([mM])(d+)[F](d+)[mM]")

def current_work():
    todo_list = [
        "NEED TO FIND JUNCTION POINT!!",
        "Strategy: only work with primary line",
        "Make synthetic data to check current code!!!",
        "Need to handle indels better (max indel size?)"
        ]
    todo_list.insert(0, "="*80)
    todo_list.append("="*80)
    
    for todo in todo_list:
        print >> sys.stderr, "*", todo
        print "*", todo

def main():
    parser = argparse.ArgumentParser(description="Extract all reads in SAM_FILE that map to REF_NAME, or have pair that maps to it (including fusion matches)")
    parser.add_argument("SAM_FILE", type=file, help="")
    # parser.add_argument("REF_NAME", help="Extract all reads w")
    parser.add_argument("--junction", metavar="JUNCTION_FILE",
                        type=argparse.FileType('w'),
                        help="Output junction table to %(metavar)s")
    parser.add_argument("--fusionreads", metavar="READ_FILE",
                        type=argparse.FileType('w'),
                        help="Output reads containing junctions to %(metavar)s")
    # parser.add_argument("--outdir", metavar="OUTDIR", help="Save output file to %(metavar)s (default is same directory as SAM_FILE)")
    # parser.add_argument("--verbose", help="increase output verbosity",action="store_true",default=False)
    # parser.add_argument("--group", help="group reads and output groups to separate FASTQ files.",action="store_true",default=False)
    args = parser.parse_args()

    junction_dict = find_chimeric_reads(args.SAM_FILE.name)
    # find_chimeric_reads_htseq(args.SAM_FILE.name)
    if args.junction:
        for key in sorted(junction_dict.keys()):
            (l_chrom, l_junc),(r_chrom, r_junc) = key
            value = junction_dict[key]
            print >>args.junction, "{4}\t{0}:{1}~{2}:{3}".format(l_chrom, l_junc,
                                                                 r_chrom, r_junc,
                                                                 len(value))
    if args.fusionreads:
        for key in sorted(junction_dict.keys()):
            (l_chrom, l_junc),(r_chrom, r_junc) = key
            for readname in junction_dict[key]:
                print >>args.fusionreads, "{4}\t{0}:{1}~{2}:{3}".format(l_chrom, l_junc,
                                                                        r_chrom, r_junc,
                                                                        readname)

def find_chimeric_reads(sam_filename):
    junction_dict = {}
    samfile = pysam.Samfile(sam_filename)
    refs = samfile.references
    chroms_re = re.compile("("+'|'.join(samfile.references)+")-("+'|'.join(samfile.references)+")")
    
    # for aread in samfile.fetch():
    print >>sys.stderr, "aread.qname, samfile.getrname(aread.tid), aread.cigarstring, aread.aend, aread.alen, aread.pos, aread.qend, aread.qlen, aread.qstart, tags['SA']"    
    for aread in samfile:
        tags = dict(aread.tags)
        if "SA" in tags:
            sa_string = tags["SA"]
            sa_list = (sa_string[:-1] if sa_string.endswith(';') else sa_string).split(';')
            # print "-"*60, "\n<SA MATCH>", aread
            if not aread.is_paired:
                readnum = 0
            elif aread.is_read1:
                readnum = 1
            elif aread.is_read2:
                readnum = 2
            else:
                raise StandardError, "Problem with read number"
            if (aread.flag & 0x900) == 0: # this the "primary line", SO HERE
                line_type = 'prime' #primary
            elif (aread.flag & 0x800) == 0x800: 
                line_type = 'suppl' #supplementary
            elif (aread.flag & 0x100) == 0x100: 
                line_type = 'secon' #secondary
            else:
                line_type = '?'

            rname = samfile.getrname(aread.tid)
            if (aread.flag & 0x900) == 0 : #read is primary
                print "{0.qname:20} r{3} {4}\t{1}:{0.pos}-{0.aend}[{0.alen}]\tQ:{0.qstart}-{0.qend}[{0.qlen}]\t{0.cigarstring:10}\tSA<{2}>".format(aread,rname, sa_list,readnum,line_type)
                if aread.is_reverse:
                    strand = "-"
                else:
                    strand = "+"
                read_parts = [ReadFragment(aread.qstart,aread.qend,aread.pos,aread.aend,rname,strand)]
                for supline in sa_list:
                    print supline
                    rname,pos,strand,CIGAR,mapQ,NM = supline.split(",")
                    pos_b0 = int(pos)-1
                    if CIGAR.count("M") != 1:
                        print "Number of matches != 1 " + CIGAR
                    for op in HTSeq.parse_cigar(CIGAR, pos_b0, rname, strand): 
                        print op, op.query_from, op.query_to, op.ref_iv
                        if op.type == "M":
                            read_parts.append(ReadFragment(op.query_from, op.query_to,op.ref_iv.start,op.ref_iv.end,op.ref_iv.chrom,strand))
                read_parts.sort()
                for part in read_parts:
                    print part
                for i in range(len(read_parts)-1):
                    l_part = read_parts[i]
                    r_part = read_parts[i+1]
                    l_chrom = l_part.rname
                    r_chrom = r_part.rname
                    if (l_part.strand == "-") and (r_part.strand == "-"):
                        l_junc = l_part.start_d
                        r_junc = r_part.end_d
                    else:
                        l_junc = l_part.end_d
                        r_junc = r_part.start_d
                    print "junction: {0.rname}:{2}<->{1.rname}:{3}".format(l_part,r_part,l_junc,r_junc)
                    junction_dict.setdefault(((l_chrom, l_junc),(r_chrom, r_junc)),set()).add(aread.qname)
                print ""
            else:
                print "{0.qname:20} r{3} {4}\t{1}:{0.pos}-{0.aend}[{0.alen}]\tQ:{0.qstart}-{0.qend}[{0.qlen}]\t{0.cigarstring:10}".format(aread,rname, sa_list,readnum,line_type)
    return junction_dict


def parse_samfile_htseq(sam_filename,refname):
    sambase,samext = os.path.splitext(sam_filename)
    if samext == ".sam":
        align_seq = iter(HTSeq.SAM_Reader( sam_filename ))
    elif samext == ".bam":
        align_seq = iter(HTSeq.BAM_Reader( sam_filename ))
    else:
        print >>sys.stderr, "Problem with SAM/BAM File:", sam_filename
        sys.exit(1)

    record_list = []
    for aread in align_seq:
        if aread.iv.chrom == refname:
            print "READ:", aread
        else:
            optional_fields = dict(aread.optional_fields)
            print aread.mate_start, optional_fields.get("XP","no xp"), optional_fields.get("XF","no xf"), aread.read.name
    
    return record_list

class ChimeraJunction:
    def __init__(self,l_part,r_part,qname):
        self.qname = qname
        l_chrom = l_part.rname
        r_chrom = r_part.rname
        if (l_part.strand == "-") and (r_part.strand == "-"):
            l_junc = l_part.start_d
            l_end = l_part.end_d
            r_junc = r_part.end_d
            r_end = r_part.start_d
        else:
            l_junc = l_part.end_d
            r_junc = r_part.start_d

            l_end = l_part.start_d
            r_end = r_part.end_d

        if l_part.rname < r_part.rname:
            self.first_chrom, self.first_junc, self.first_end = l_chrom, l_junc, l_end
            self.second_chrom, self.second_junc, self.second_end = r_chrom, r_junc, r_end
        else:
            self.first_chrom, self.first_junc, self.first_end = r_chrom, r_junc, r_end
            self.second_chrom, self.second_junc, self.second_end = l_chrom, l_junc, l_end
           
        print "junction: {0.rname}:{2}<->{1.rname}:{3}".format(l_part,r_part,l_junc,r_junc)
        junction_dict.setdefault(((l_chrom, l_junc),(r_chrom, r_junc)),set()).add(aread.qname)

    def __str__(self):
        return "{0}:{1}~{2}:{3}".format(first_chrom, first_junc,second_chrom, second_junc)

    def __lt__(self,other):
        if self.first_chrom < other.first_chrom:
            return True
        elif self.first_chrom > other.first_chrom:
            return False
        elif self.first_chrom == other.first_chrom:
            if self.second_chrom < other.second_chrom:
                return True
            elif self.second_chrom > other.second_chrom:
                return False
            elif self.second_chrom == other.second_chrom:
                if self.first_junc < other.first_junc:
                    return True
                elif self.first_junc > other.first_junc:
                    return False
                elif self.first_junc == other.first_junc:
                    if self.second_junc < other.second_junc:
                        return True
                    elif self.second_junc > other.second_junc:
                        return False
                    elif self.second_junc == other.second_junc:
                        if self.first_end < other.first_end:
                            return True
                        elif self.first_end > other.first_end:
                            return False
                        elif self.first_end == other.first_end:
                            if self.second_end < other.second_end:
                                return True
                            else:
                                return False
        else:
            raise StandardError, "Problem in chrom comparison"



class ReadFragment:
    def __init__(self,qstart,qend,pos,aend,rname,strand):
        self.qstart = qstart
        self.qend = qend
        self.pos = pos+1
        self.aend = aend
        self.rname = rname
        self.strand = strand
        if strand == "-":
            self.start_d = self.aend
            self.end_d = self.pos
        else:
            self.start_d = self.pos
            self.end_d = self.aend
    def __lt__(self,other):
        return self.qstart < other.qstart

    def __str__(self):
        return "Q:{0.qstart}-{0.qend} -> {0.rname}:{0.pos}-{0.aend}{0.strand} ({0.start_d}-{0.end_d})".format(self)
    
    # @property
    # def range(self):
    #     return "{0.chrom}:{0.start}-{0.end}".format(self.iv)
    #     # return self.iv

    # @property
    # def count(self):
    #     return len(self.read_list)

    # @property
    # def insertion_point(self):
    #     if self._insert_point == None:
    #         self._determine_insertion_point()
    #     return self._insert_point

    # @property
    # def insertion_side(self):
    #     if self._insert_side == None:
    #         self._determine_insertion_point()
    #     return self._insert_side


if __name__ == "__main__":
    current_work()
    main()
    current_work()

