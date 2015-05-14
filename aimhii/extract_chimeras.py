#!/usr/bin/env python

import HTSeq
import argparse
import sys
import pysam
import re
import os
from Bio import SeqIO
from find_insert_clusters import ReadCluster, filter_clusters, find_cluster_pairs, RIGHT, LEFT, ClusterGroup
import csv
from operator import attrgetter

from plot_metaclusters import plot_clusters, READ_PITCH, READ_HEIGHT
import matplotlib

import logging
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

from aimhii.__about__ import (
    __author__, __copyright__, __email__, __license__, __summary__, __title__,
    __uri__, __version__)

# cigar_re = re.compile("(d+)([mM])(d+)[F](d+)[mM]")
OPPOSITE_STRAND = {"+":"-", "-":"+"}
OPPOSITE_SIDE = {RIGHT:LEFT, LEFT:RIGHT}
MINUS = "-";PLUS = "+"

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
    parser.add_argument("--table", metavar="TABLE_FILE", default=sys.stdout,
                        type=argparse.FileType('w'),
                        help="Output detailed table of insertions to %(metavar)s")
    # parser.add_argument("--cluster", action="store_true", help="Cluster primary_fragments of junction reads", default=False)
    # parser.add_argument("--insertonly", action="store_true", help="Filter out junctions that don't include the insert", default=False)
    parser.add_argument("--alljunctions", action="store_true", help="Process all junctions (default: if --insert is provide, exclude reads that don't contain the insert", default=False)
    parser.add_argument("--minreads", type=int, metavar="NUM_READS", default=2, 
                        help="Clusters consisting of less than %(metavar)s reads are ignored (default: %(default)s)")
    parser.add_argument("--maxgap", type=int, metavar="GAP_LENGTH", default=5000, 
                        help="Clusters separated by no more than %(metavar)s) bases are considered to be cluster doublets (default: %(default)s)")
    parser.add_argument("--plot", metavar="PLOT_FILE_PREFIX",
                        help="Generate plots of metaclusters.  Plot file names will be prefixed with %(metavar)s")
    parser.add_argument('-V', '--version', action='version', version="%(prog)s v{0}".format(__version__))

    # parser.add_argument("--outdir", metavar="OUTDIR", help="Save output file to %(metavar)s (default is same directory as SAM_FILE)")
    # parser.add_argument("--verbose", help="increase output verbosity",action="store_true",default=False)
    # parser.add_argument("--group", help="group reads and output groups to separate FASTQ files.",action="store_true",default=False)
    args = parser.parse_args()

    junction_list, metacluster_list = run_analysis(args.SAM_FILE.name, args.minreads,
                                                   args.maxgap, args.table,
                                                   args.insert, not(args.alljunctions))
    #--------------------------------------------------
    if args.junction:
        output_junctions(junction_list, args.junction)
        
    if args.fusionreads:
        output_fusionreads(junction_list, args.fusionreads)

    if args.plot:
        plot_clusters(metacluster_list,args.plot)
    #--------------------------------------------------
    
def output_junctions(junction_list,junction_handle):
    # count the number of hits at each junction
    junction_count_dict = {}
    for curjunc in junction_list:
        junction_count_dict.setdefault(curjunc.junction_tuple,[]).append(curjunc)
        # junction_count_dict[curjunc.junction_tuple] = 1+junction_count_dict.get(curjunc.junction_tuple,0)
    for key in sorted(junction_count_dict.keys()):
        l_chrom, l_junc,r_chrom, r_junc = key
        print >>junction_handle, "{4}\t{0}:{1}~{2}:{3}".format(l_chrom, l_junc,
                                                             r_chrom, r_junc,
                                                             len(junction_count_dict[key]))

def output_fusionreads(junction_list, fusion_handle):
    # for curjunc in sorted(junction_list,key=lambda x: x.junction_tuple):
    for curjunc in junction_list:
        print >>fusion_handle, "{0.readname}\t{0}".format(curjunc)

def run_analysis(mapread_bam, minreads, maxgap, outhandle, insert_seq=None, insertonly=False):
    if insert_seq:
        insert_rec = SeqIO.read(insert_seq,"fasta")
        ChimeraJunction.InsertSeqID = insert_rec.id

    junction_list = find_chimeric_reads(mapread_bam)
    if insertonly and insert_seq:
        junction_list = [junction for junction in junction_list if junction.contains_insert]
    junction_list.sort()
    # find_chimeric_reads_htseq(args.SAM_FILE.name)

    
    
    logger.debug("-"*20 + "CLUSTERING" + "-"*20)
    cluster_list = cluster_junctions(junction_list)
    cluster_list = filter_clusters(cluster_list, minreads)
    cluster_list.sort()

    for cluster in cluster_list:
        logger.debug(cluster)

    for cluster in cluster_list:
        logger.debug("{0.count:<10} {0.range}".format(cluster))

    for cluster in cluster_list:
        logger.debug("DETAILS {0}".format(cluster.insert_details))


    metacluster_list = find_cluster_pairs(cluster_list,maxgap)
    for metacluster in metacluster_list:
        logger.debug(metacluster.str_with_secondary)

    outwriter = csv.writer(outhandle,dialect=csv.excel)
    # print >>sys.stderr, "metacluster_list", metacluster_list
    outwriter.writerow(ClusterGroup.Header)
    for metacluster in metacluster_list:
        # print metacluster.output
        outwriter.writerow(metacluster.output)

    return junction_list, metacluster_list

def find_chimeric_reads(sam_filename):
    junction_list = []
    samfile = pysam.Samfile(sam_filename)
    # refs = samfile.references
    # chroms_re = re.compile("("+'|'.join(samfile.references)+")-("+'|'.join(samfile.references)+")")
    
    # for aread in samfile.fetch():
    logger.debug("aread.qname, samfile.getrname(aread.tid), aread.cigarstring, aread.aend, aread.alen, aread.pos, aread.qend, aread.qlen, aread.qstart, tags['SA']")
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
                if aread.is_reverse:
                    primary_strand = MINUS
                else:
                    primary_strand = PLUS
                logger.debug("{0.qname:20} r{3} {4}\t{1}:{0.pos}-{0.aend}[{0.alen}]\tQ:{0.qstart}-{0.qend}[{0.qlen}]\t{5}{0.cigarstring:10}\tSA<{2}>".format(aread,rname, sa_list,readnum,line_type,primary_strand))
                read_parts = ReadFragment.listFromCIGAR(aread.cigarstring, aread.pos, rname, primary_strand)
                for supline in sa_list:
                    logger.debug(supline)
                    rname,pos,sup_strand,CIGAR,mapQ,NM = supline.split(",")
                    pos_b0 = int(pos)-1
                    # if CIGAR.count("M") != 1:
                    #     print "Number of matches != 1 ", CIGAR, "SKIPPING THIS ONE!!!", "LOOK INTO MERGING ACROSS DELETES?"
                    #     continue
                    read_parts.extend(ReadFragment.listFromCIGAR(CIGAR, pos_b0, rname, sup_strand))
                read_parts.sort(key=attrgetter('qstart'))
                for part in read_parts:
                    logger.debug(part)
                for i in range(len(read_parts)-1):
                    cur_junction = ChimeraJunction(read_parts[i],read_parts[i+1],aread.qname)
                    logger.debug("{0} {1}".format(cur_junction, cur_junction.readname))
                    junction_list.append(cur_junction)
                logger.debug("."*10)
            else:
                logger.debug("{0.qname:20} r{3} {4}\t{1}:{0.pos}-{0.aend}[{0.alen}]\tQ:{0.qstart}-{0.qend}[{0.qlen}]\t{0.cigarstring:10}".format(aread,rname, sa_list,readnum,line_type))
    return junction_list

def cluster_junctions(junction_list):
    cluster_list = []
    for cur_junction in junction_list:
        added_to_cluster = False
        for cur_cluster in cluster_list:
            if cur_cluster.matches_cluster(cur_junction):
                cur_cluster.add_read(cur_junction)
                added_to_cluster = True
        if not added_to_cluster:
            # aread doesn't overlap cur_cluster, so time to start a new cluster
            cluster_list.append(ReadCluster(cur_junction))
    return cluster_list

class ChimeraJunction:
    InsertSeqID = None
    def __init__(self,l_part,r_part,qname):
        self.readname = qname

        ##========================
        ## Select primary fragment
        ##========================
        if (((l_part.rname == ChimeraJunction.InsertSeqID) or
            (r_part.rname == ChimeraJunction.InsertSeqID)) and
            (l_part.rname != r_part.rname)):
            # Organize relative to reference chrom, not insert
            if l_part.rname == ChimeraJunction.InsertSeqID:
                self.primary_frag = r_part # ref_part = r_part
                self.secondary_frag = l_part # insert_part = l_part
            else:
                self.primary_frag = l_part # insert_part
                self.secondary_frag = r_part # ref_part
        elif l_part.rname < r_part.rname or (l_part.rname == r_part.rname and
                                             (l_part.start < r_part.start or
                                              (l_part.start == r_part.start and
                                               l_part.end < r_part.end))):
            self.primary_frag = l_part
            self.secondary_frag = r_part
        elif l_part.rname > r_part.rname or (l_part.rname == r_part.rname and
                                             (l_part.start > r_part.start or
                                             (l_part.start == r_part.start and
                                               l_part.end > r_part.end))):
            self.primary_frag = r_part
            self.secondary_frag = l_part
        else:
            print "STRANGE READ!!! NEED TO CHECK!!! l_part.rname == r_part.rname and l_part.start == r_part.start??"
            print "qname:{0}\tl_part:{1}\tr_part:{2}".format(qname,l_part,r_part)
            self.primary_frag = l_part
            self.secondary_frag = r_part

        ##========================
        ## Normalize Junctions:
        ## primary_frag should always be on plus strand, so invert junctions where primary_frag is minus strand
        ##========================
        if self.primary_frag.strand == PLUS and self.secondary_frag.strand == PLUS :
            l_part.distal = l_part.start+1
            l_part.junc = l_part.end
            r_part.junc = r_part.start+1
            r_part.distal = r_part.end
        elif self.primary_frag.strand == MINUS and self.secondary_frag.strand == MINUS :
            self.primary_frag.strand = self.secondary_frag.strand = PLUS
            l_part,r_part = r_part,l_part
            l_part.distal = l_part.start+1
            l_part.junc = l_part.end
            r_part.junc = r_part.start+1
            r_part.distal = r_part.end
        elif self.primary_frag.strand == PLUS and self.secondary_frag.strand == MINUS :
            if self.primary_frag == l_part:
                l_part.distal = l_part.start+1
                l_part.junc = l_part.end
                r_part.junc = r_part.end
                r_part.distal = r_part.start+1
            else: # self.primary_frag == r_part:
                l_part.distal = l_part.end
                l_part.junc = l_part.start+1
                r_part.junc = r_part.start+1
                r_part.distal = r_part.end
        elif self.primary_frag.strand == MINUS and self.secondary_frag.strand == PLUS :
            self.primary_frag.strand,self.secondary_frag.strand = PLUS,MINUS
            l_part,r_part = r_part,l_part
            if self.primary_frag == l_part:
                l_part.distal = l_part.start+1
                l_part.junc = l_part.end
                r_part.junc = r_part.end
                r_part.distal = r_part.start+1
            else: # self.primary_frag == r_part:
                l_part.distal = l_part.end
                l_part.junc = l_part.start+1
                r_part.junc = r_part.start+1
                r_part.distal = r_part.end
        else:
            raise StandardError, "Problem with strand combinations"

        self.l_frag = l_part
        self.r_frag = r_part
        if self.primary_frag == l_part:
            self.insert_side = RIGHT
        else:
            self.insert_side = LEFT
            
        self.insert_point = self.l_frag.junc

        logger.debug("l.distal:{0.distal} l.junc:{0.junc} r.junc:{1.junc} r.distal:{1.distal}".format(l_part,r_part))
        logger.debug("l.start:{0.start} l.end:{0.end} r.start:{1.start} r.end:{1.end}".format(l_part,r_part))

        self.iv = self.primary_frag
           
        logger.debug("junction: {0.primary_frag.rname}:{0.primary_frag.junc}<->{0.secondary_frag.rname}:{0.secondary_frag.junc} insert_side:{0.insert_side}".format(self))

    def __str__(self):
        return "{0.l_frag.chrom}:{0.l_frag.junc}~{0.r_frag.chrom}:{0.r_frag.junc} insert:{0.insert_side}".format(self)

    # def __hash__(self):
    #     # return hash((self.l_chrom, self.l_junc),(self.r_chrom, self.r_junc))
    #     return hash(((self.chrom1, self.junc1),(self.chrom2, self.junc2)))
    @property
    def junction_tuple(self):
        return (self.l_frag.chrom, self.l_frag.junc, self.r_frag.chrom, self.r_frag.junc)

    @property
    def contains_insert(self):
        return (ChimeraJunction.InsertSeqID == self.primary_frag.rname or
                ChimeraJunction.InsertSeqID == self.secondary_frag.rname)

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
        elif self.primary_frag.distal < other.primary_frag.distal:
            return True
        elif self.primary_frag.distal > other.primary_frag.distal:
            return False
        elif self.secondary_frag.distal < other.secondary_frag.distal:
            return True
        elif self.secondary_frag.distal > other.secondary_frag.distal:
            return False
        else:
            return self.readname < other.readname
    def draw(self,i):
        xy = (self.primary_frag.start,i*READ_PITCH)
        width = self.primary_frag.length
        height = READ_HEIGHT
        logger.debug("I am pretending to draw myself xy:{0}, width:{1}, height:{2}".format(xy,width,height))
        if self.insert_side == RIGHT:
            rect_color = "green"
        else:
            rect_color = "purple"        
        rect = matplotlib.patches.Rectangle( xy, width=width, height=height,facecolor=rect_color,edgecolor="none",alpha=0.6)
        # scale_factor = 0.4
        # if abs(dx) < scaled_headlen:
        #     headlen = 0
        # else:
        #     headlen = scaled_headlen
        # exon_arrow = ax.arrow(x1,y,dx,0, shape='full', width=exon_dy, length_includes_head=True, head_width=scale_factor,head_length=headlen, transform=trans,fc="black",alpha=0.9)

        return rect

class ReadFragment(HTSeq.GenomicInterval):
    def __init__(self,qstart,qend,pos,aend,rname,strand):
        self.qstart = qstart
        self.qend = qend
        self.rname = rname
        super(ReadFragment,self).__init__(rname,pos,aend,strand)
        logger.debug("SUPER", super(ReadFragment,self).__str__())

    def __lt__(self,other):
        # this is important for ordering fragments in find_chimeric_reads()
        # return self.qstart < other.qstart
        raise NotImplementedError

    def __str__(self):
        return "Q:{0.qstart}-{0.qend} -> {0.rname}:{0.start}-{0.end}{0.strand} ({0.start_d}-{0.end_d})".format(self)

    def extend(self,qstart,qend,pos,aend,rname,strand):
        """Extend ReadFragment to include portion beyond a delete"""
        if self.rname != rname:
            raise StandardError, "Incompatible reference sequences: {0} != {1}".format(self.rname, rname)
        elif self.strand != strand:
            raise StandardError, "Incompatible strands: {0} != {1}".format(self.strand, strand)
        else:
            logger.debug("Before extension: {0}".format(self))
            self.qstart = min(self.qstart,qstart)
            self.qend = max(self.qend,qend)
            new_iv = HTSeq.GenomicInterval(rname,pos,aend,strand)
            logger.debug("new_iv:", new_iv)
            self.extend_to_include(new_iv)
            logger.debug("After extension: {0}".format(self))
        # raise NotImplementedError

    @classmethod
    def listFromCIGAR(cls, cigarstring,position_b0, refname, strand):
        read_parts = []
        if strand == MINUS: # need to reverse the CIGAR
            logger.debug("Reversing CIGAR for minus strand read fragment")
            cigarstring = "".join(reversed(re.findall("\d+[MIDNSHP=X]", cigarstring)))

        op_type_list = []
        for op in HTSeq.parse_cigar(cigarstring, position_b0, refname, strand):
            logger.debug(map(str,(op, op.query_from, op.query_to, op.ref_iv)))
            if op.type == "M":
                if "M" in op_type_list:
                    if len(op_type_list) >=2 and op_type_list[-1] == "D" and op_type_list[-2] == "M":
                        logger.debug(map(str,("extending (D):", op, op.query_from, op.query_to, op.ref_iv)))
                        read_parts[-1].extend(op.query_from, op.query_to,op.ref_iv.start,op.ref_iv.end,op.ref_iv.chrom,strand)
                    elif len(op_type_list) >=2 and op_type_list[-1] == "I" and op_type_list[-2] == "M":
                        logger.debug(map(str,("extending (I):", op, op.query_from, op.query_to, op.ref_iv)))
                        read_parts[-1].extend(op.query_from, op.query_to,op.ref_iv.start,op.ref_iv.end,op.ref_iv.chrom,strand)
                    else:
                        logger.debug("CIGAR WARNING: Number of matches > 1: {0}".format(cigarstring))
                else:
                    logger.debug(map(str,("appending:", op, op.query_from, op.query_to, op.ref_iv)))
                    suppl_frag = cls(op.query_from, op.query_to,op.ref_iv.start,op.ref_iv.end,op.ref_iv.chrom,strand)
                    read_parts.append(suppl_frag)
            op_type_list.append(op.type)
        return read_parts

if __name__ == "__main__":
    main()
    # logger.debug('debug message')
    # logger.info('info message')
    # logger.warn('warn message')
    # logger.error('error message')
    # logger.critical('critical message')
