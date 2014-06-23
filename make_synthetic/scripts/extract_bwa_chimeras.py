import HTSeq
import argparse
import sys
import pysam
import re
import os
from Bio import SeqIO

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
    parser.add_argument("SAM_FILE", type=file, help="BAM/SAM file containing mapped reads")
    parser.add_argument("--insert", type=file, help="The sequence of the inserted DNA fragment.  Used to determine how to organize junctions.")
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

    if args.insert:
        insert_rec = SeqIO.read(args.insert,"fasta")
        ChimeraJunction.InsertSeqID = insert_rec.id

    junction_list = find_chimeric_reads(args.SAM_FILE.name)
    # find_chimeric_reads_htseq(args.SAM_FILE.name)
    if args.junction:
        # count the number of hits at each junction
        junction_count_dict = {}
        for curjunc in junction_list:
            junction_count_dict[curjunc.junction_tuple] = 1+junction_count_dict.get(curjunc.junction_tuple,0)
        for key in sorted(junction_count_dict):
            (l_chrom, l_junc),(r_chrom, r_junc) = key
            print >>args.junction, "{4}\t{0}:{1}~{2}:{3}".format(l_chrom, l_junc,
                                                                 r_chrom, r_junc,
                                                                 junction_count_dict[key])
    if args.fusionreads:
        # for curjunc in sorted(junction_list,key=lambda x: x.junction_tuple):
        for curjunc in sorted(junction_list):
            print >>args.fusionreads, "{0.readname}\t{0.chrom1}:{0.junc1}~{0.chrom2}:{0.junc2}".format(curjunc)

def find_chimeric_reads(sam_filename):
    junction_list = []
    samfile = pysam.Samfile(sam_filename)
    # refs = samfile.references
    # chroms_re = re.compile("("+'|'.join(samfile.references)+")-("+'|'.join(samfile.references)+")")
    
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
                    cur_junction = ChimeraJunction(read_parts[i],read_parts[i+1],aread.qname)
                    # l_chrom = l_part.rname
                    # r_chrom = r_part.rname
                    # if (l_part.strand == "-") and (r_part.strand == "-"):
                    #     l_junc = l_part.start_d
                    #     r_junc = r_part.end_d
                    # else:
                    #     l_junc = l_part.end_d
                    #     r_junc = r_part.start_d
                    # print "junction: {0.rname}:{2}<->{1.rname}:{3}".format(l_part,r_part,l_junc,r_junc)
                    # junction_dict.setdefault(((l_chrom, l_junc),(r_chrom, r_junc)),set()).add(aread.qname)
                    print cur_junction
                    # junction_dict.setdefault(cur_junction,set()).add(aread.qname)
                    junction_list.append(cur_junction)
                print ""
            else:
                print "{0.qname:20} r{3} {4}\t{1}:{0.pos}-{0.aend}[{0.alen}]\tQ:{0.qstart}-{0.qend}[{0.qlen}]\t{0.cigarstring:10}".format(aread,rname, sa_list,readnum,line_type)
    return junction_list


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
    InsertSeqID = None
    def __init__(self,l_part,r_part,qname):
        self.readname = qname
        # self.l_chrom = l_part.rname
        # self.r_chrom = r_part.rname
        l_chrom = l_part.rname
        r_chrom = r_part.rname
        if (l_part.strand == "-") and (r_part.strand == "-"):
            l_part.junc = l_part.start_d
            l_part.end = l_part.end_d
            r_part.junc = r_part.end_d
            r_part.end = r_part.start_d
        else:
            l_part.junc = l_part.end_d
            r_part.junc = r_part.start_d
            l_part.end = l_part.start_d
            r_part.end = r_part.end_d

        if (l_part.rname == ChimeraJunction.InsertSeqID) or (r_part.rname == ChimeraJunction.InsertSeqID) and (l_part.rname != r_part.rname):
            # Organize in the relative to reference chrom, not insert
            if l_part.rname == ChimeraJunction.InsertSeqID:
                self.primary_frag = r_part # ref_part = r_part
                self.secondary_frag = l_part # insert_part = l_part
            else:
                self.primary_frag = l_part # insert_part
                self.secondary_frag = r_part # ref_part
            # self.primary_frag = ref_part
            # self.secondary_frag = insert_part
        elif l_part.rname < r_part.rname or (l_part.rname == r_part.rname and l_part < r_part):
            self.primary_frag = l_part
            self.secondary_frag = r_part
        elif l_part.rname > r_part.rname or (l_part.rname == r_part.rname and l_part > r_part):
            self.primary_frag = r_part
            self.secondary_frag = l_part
        else:
            raise StandardError
        
        if self.primary_frag.junc > self.primary_frag.end:
            # junction is at right end of ref fragment
            self.chrom1, self.junc1, self.end1 = self.primary_frag.rname, self.primary_frag.junc, self.primary_frag.end
            self.chrom2, self.junc2, self.end2 = self.secondary_frag.rname, self.secondary_frag.junc, self.secondary_frag.end
        else:
            # junction is at left end of ref fragment
            self.chrom1, self.junc1, self.end1 = self.secondary_frag.rname, self.secondary_frag.junc, self.secondary_frag.end
            self.chrom2, self.junc2, self.end2 = self.primary_frag.rname, self.primary_frag.junc, self.primary_frag.end
        # elif l_part.rname < r_part.rname:
        #     self.chrom1, self.junc1, self.end1 = l_chrom, l_part.junc, l_part.end
        #     self.chrom2, self.junc2, self.end2 = r_chrom, r_part.junc, r_part.end
        # else:
        #     self.chrom1, self.junc1, self.end1 = r_chrom, r_part.junc, r_part.end
        #     self.chrom2, self.junc2, self.end2 = l_chrom, l_part.junc, l_part.end
           
        print "junction: {0.rname}:{0.junc}<->{1.rname}:{1.junc}".format(l_part,r_part)
        # junction_dict.setdefault(((l_chrom, l_junc),(r_chrom, r_junc)),set()).add(aread.qname)

    def __str__(self):
        return "{0.chrom1}:{0.junc1}~{0.chrom2}:{0.junc2}".format(self)

    # def __hash__(self):
    #     # return hash((self.l_chrom, self.l_junc),(self.r_chrom, self.r_junc))
    #     return hash(((self.chrom1, self.junc1),(self.chrom2, self.junc2)))
    @property
    def junction_tuple(self):
        return (self.chrom1, self.junc1),(self.chrom2, self.junc2)

    def __lt__(self,other):
        if self.primary_frag.rname < other.primary_frag.rname:
            return True
        elif self.primary_frag.rname > other.primary_frag.rname:
            return False
        elif self.secondary_frag.rname < other.secondary_frag.rname:
            return True
        elif self.secondary_frag.rname > other.secondary_frag.rname:
            return False
        #--------------------------------------------------
        elif self.primary_frag.junc < other.primary_frag.junc:
            return True
        elif self.primary_frag.junc > other.primary_frag.junc:
            return False
        elif self.secondary_frag.junc < other.secondary_frag.junc:
            return True
        elif self.secondary_frag.junc > other.secondary_frag.junc:
            return False
        #--------------------------------------------------
        elif self.primary_frag.end < other.primary_frag.end:
            return True
        elif self.primary_frag.end > other.primary_frag.end:
            return False
        elif self.secondary_frag.end < other.secondary_frag.end:
            return True
        elif self.secondary_frag.end > other.secondary_frag.end:
            return False
        else:
            return self.readname < other.readname


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

