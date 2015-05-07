from __future__ import division
import HTSeq
import argparse
import sys
# import re
import os
import numpy as np
# import csv
# import collections
# import itertools
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio import SeqIO
# ambiguous_re = re.compile("ambiguous\[(.*)\]")
# sample_re = re.compile("(.*)_[ACGT]{6}_L002_R1_001")

import logging
logger = logging.getLogger(__name__)
# logger.setLevel(logging.DEBUG)
logger.setLevel(logging.INFO)

LEFT="left"
RIGHT="right"
CLUSTER_GFF_TYPE="experimental_result_region"
GAP_GFF_TYPE="insertion_site"
SINGLETON_INSERT_GFF_TYPE="junction"
"""
python2.7 $SCRIPTS/find_insert_clusters.py $COLLAB/AlspaughLab/clipped_mapped/PE_H99align_VSL_093013.clips.sort.bam > $COLLAB/AlspaughLab/clipped_mapped/clusters_v1.txt
"""

def main():
    parser = argparse.ArgumentParser(description="Identifies read clusters (sets of overlapping reads), MORE")
    parser.add_argument("SAM_FILE", type=file, help="")
    parser.add_argument("--maxgap", type=int, metavar="GAP_LENGTH", default=100, 
                        help="Clusters separated by no more than %(metavar)s) bases are considered to be cluster doublets (default: %(default)s)")
    parser.add_argument("--minreads", type=int, metavar="NUM_READS", default=2, 
                        help="Clusters consisting of less than %(metavar)s reads are ignored (default: %(default)s)")
    parser.add_argument("--gff",type=argparse.FileType('w'),metavar="GFF_FILE",
                        help="Output results in GFF format to %(metavar)s")
    # parser.add_argument("--outdir", metavar="OUTDIR", help="Save output file to %(metavar)s (default is same directory as SAM_FILE)")
    # parser.add_argument("--verbose", help="increase output verbosity",action="store_true",default=False)
    args = parser.parse_args()

    
    # out_fastq_name = os.path.splitext(args.SAM_FILE.name)[0] + ".clips.fastq"
    # if args.outdir:
    #     out_fastq_name = os.path.join(args.outdir,os.path.basename(out_fastq_name))
    print "NEEED TO MAKE VALUES base 1!!!!", "#"*50
    read_iter = load_sam_or_bam(args.SAM_FILE.name)
    cluster_list = find_clusters(read_iter)
    cluster_list = filter_clusters(cluster_list, args.minreads)
    
    for cluster in cluster_list:
        print "{0.count:<10} {0.range}".format(cluster)

    print "NEEED TO MAKE VALUES base 1!!!!", "#"*50
    cluster_group_list = find_cluster_pairs(cluster_list,args.maxgap)
    # SeqIO.write(record_list, out_fastq_name, "fastq-sanger")
    print "NEEED TO MAKE VALUES base 1!!!!", "#"*50
    for group in cluster_group_list:
        print group

    if args.gff:
        generate_gff(args.gff,cluster_group_list)
    # print >>sys.stderr, "="*60
    # for match in match_list:
    #     print match
    #     # print >>sys.stderr, "MATCH:", match.get_gff_line()
    print "NEEED TO MAKE VALUES base 1!!!!", "#"*50


def generate_gff(gff_outhandle, cluster_group_list):
    print "NEED TO FIGURE OUT INSERT STRAND"
    for i, group in enumerate(cluster_group_list):
        for feature in group.get_features(i+1):
            gff_outhandle.write(feature.get_gff_line())
    """
            /Users/josh/Documents/BioinfCollabs/docs/sofa.obo
    name: region
name: read_pair
name: match_set
name: match_part
name: read
name: flanking_region
name: clip
name: tag
name: match
name: insertion_site
name: sequence_difference
name: junction
name: remark
name: experimental_result_region
name: gap
name: sequence_alteration
name: biomaterial_region
name: experimental_feature
name: paired_end_fragment
name: gene_group
name: adjacent_to
name: associated_with
name: connects_on
name: contained_by
name: contains
name: derives_from
name: disconnected_from
name: edited_from
name: edited_to
name: lost
name: overlaps
name: position_of
            """
def filter_clusters(cluster_list, minreads):
    return [cluster for cluster in cluster_list if cluster.count >= minreads]
    
def find_cluster_pairs(cluster_list,max_gap_size):
    cluster_list.sort()
    metacluster_list = []
    i = 0
    while i < len(cluster_list):
        cur_cluster = cluster_list[i]
        if i == len(cluster_list)-1:
            metacluster_list.append(ClusterSingleton(cur_cluster))
            logger.debug("Last cluster is a singleton: {0}".format(metacluster_list[-1]))
            break
        else:
            next_cluster = cluster_list[i+1]

        logger.debug("cur_cluster:{0}\tnext_cluster:{1}".format(cur_cluster, next_cluster))
        if (cur_cluster.iv.chrom != next_cluster.iv.chrom):
            logger.debug("Last cluster on chrom: {0}".format(cur_cluster))
            metacluster_list.append(ClusterSingleton(cur_cluster))
            logger.debug("Found a cluster singleton: {0}".format(metacluster_list[-1]))
            i += 1
        elif cur_cluster.insertion_side == RIGHT and next_cluster.insertion_side == LEFT:
            logger.debug("Clusters face each other: {0} {1}".format(cur_cluster, next_cluster))
            # for clust in (cur_cluster, next_cluster):
            #     print >>sys.stderr, clust
            #     for read in clust.read_list:
            #         print read
            metacluster = ClusterDoublet(cur_cluster,next_cluster)
            if metacluster.gap_length <= max_gap_size:
                metacluster_list.append(metacluster)
                logger.debug("Found a cluster doublet: {0} {1}".format(metacluster.gap_length, metacluster))
                i += 2
            else:
                logger.debug("Gap is too big: {0} {1} {2}".format(cur_cluster, next_cluster, metacluster.gap_length))
                metacluster_list.append(ClusterSingleton(cur_cluster))
                logger.debug("Found a cluster singleton: {0}".format(metacluster_list[-1]))
                i += 1
        else:
            metacluster_list.append(ClusterSingleton(cur_cluster))
            logger.debug("Found a cluster singleton: {0}".format(metacluster_list[-1]))
            i += 1
    return metacluster_list
        
def load_sam_or_bam(sam_filename):
    sambase,samext = os.path.splitext(sam_filename)
    if samext == ".sam":
        align_seq = iter(HTSeq.SAM_Reader( sam_filename ))
    elif samext == ".bam":
        align_seq = iter(HTSeq.BAM_Reader( sam_filename ))
    else:
        print >>sys.stderr, "Problem with SAM/BAM File:", sam_filename
        sys.exit(1)
    return align_seq

# def find_clusters(align_seq):
#     cluster_list = []
#     cur_cluster = None
#     for aread in align_seq:
#         if not aread.aligned:
#             continue
#         elif cur_cluster == None:
#             #prime the pump
#             cur_cluster = ReadCluster(aread)
#             cluster_list.append(cur_cluster)
#         elif cur_cluster.overlaps(aread):
#             # aread overlaps cur_cluster, so it should be added
#             cur_cluster.add_read(aread)
#             # cur_cluster.add_read(aread)
#             # cur_cluster_iv.extend_to_include(aread.iv)
#         else:
#             # aread doesn't overlap cur_cluster, so time to start a new cluster
#             # print "Finished Cluster:", cur_cluster.range, cur_cluster.count
#             cur_cluster = ReadCluster(aread)
#             cluster_list.append(cur_cluster)
#     # cluster_list.append(cur_cluster)
#     # print cur_cluster.range, cur_cluster.count
#     return cluster_list

class ReadCluster:
    def __init__(self,aread):
        self.read_list = [aread]
        self.secondary_name = aread.secondary_frag.rname
        logger.debug("||".join(map(str, (aread, aread.iv, repr(aread.iv), aread.readname))))
        self.iv = aread.iv.copy()
        # self.secondary_iv = aread.secondary_frag.iv.copy()
        self.secondary_iv = aread.secondary_frag.copy()
        logger.debug("-"*60)
        logger.debug("INIT CHECK STRAND {0}.iv {0}.secondary_frag {0}.insert_side {0}.insert_point {0}.readname".format(aread))
        self.iv.strand = self.secondary_iv.strand = "."
        self._insert_point = aread.insert_point
        self._insert_side = aread.insert_side

    def add_read(self, aread):
        logger.debug("ADD CHECK STRAND {0}.iv {0}.secondary_frag {0}.insert_side {0}.insert_point {0}.readname".format(aread))
        self.iv.extend_to_include(aread.iv)
        self.secondary_iv.extend_to_include(aread.secondary_frag)
        self.read_list.append(aread)

    def matches_cluster(self,aread):
        return (self.overlaps(aread) and
                self.secondary_name == aread.secondary_frag.rname and
                self.insertion_side == aread.insert_side)

    def overlaps(self,aread):
        return self.iv.overlaps(aread.iv)

    @property
    def range(self):
        return "{0.chrom}:{0.start}-{0.end}".format(self.iv)
        # return self.iv

    @property
    def count(self):
        return len(self.read_list)

    @property
    def insertion_point(self):
        self._determine_insertion_point()
        return self._insert_point

    @property
    def insertion_side(self):
        # if self._insert_side == None:
        #     self._determine_insertion_point()
        return self._insert_side

    def _determine_insertion_point(self):
        insert_median = int(np.median([read.insert_point for read in self.read_list]))
        self._insert_point = insert_median
        # if len(self.read_list) == 1:
        #     self._insert_point = -1
        #     self._insert_side = "?"
        # start_median = int(np.median([read.iv.start for read in self.read_list]))
        # end_median = int(np.median([read.iv.end for read in self.read_list]))

        # frac_reads_start_median = sum(read.iv.start == start_median for read in self.read_list)/len(self.read_list)
        # frac_reads_end_median = sum(read.iv.end == end_median for read in self.read_list)/len(self.read_list)

        # if (frac_reads_start_median > 0.5) and (frac_reads_end_median > 0.5):
        #     print "Can't decide which end has insert point"
        # elif (frac_reads_start_median > 0.5):
        #     self._insert_point = start_median
        #     self._insert_side = LEFT
        # elif (frac_reads_end_median > 0.5):
        #     self._insert_point = end_median
        #     self._insert_side = RIGHT
        # else:
        #     self._insert_point = -2
        #     self._insert_side = "?!!"
        
        # # all_reads_start_median = all(read.iv.start == start_median for read in self.read_list)/len(self.read_list)
        # # all_reads_end_median = all(read.iv.end == end_median for read in self.read_list)/len(self.read_list)

        # # if all_reads_start_median and all_reads_end_median:
        # #     print "Can't decide which end has insert point"
        # # elif all_reads_start_median:
        # #     self._insert_point = start_median
        # #     self._insert_side = LEFT
        # # elif all_reads_end_median:
        # #     self._insert_point = end_median
        # #     self._insert_side = RIGHT
        # # else:
        # #     self._insert_point = -1
        # #     self._insert_side = "?"

    def __lt__(self,other):
        if self.iv.chrom == other.iv.chrom:
            return self.iv.start < other.iv.start
        else:
            return self.iv.chrom < other.iv.chrom

    # def distance(self,other):
    #     return max(self.iv.start-other.iv.end, other.iv.start-self.iv.end)

    def __str__(self):
        return "{0.iv.chrom}:{0.iv.start}-{0.iv.end} LENGTH:{0.iv.length} COUNT:{0.count} INSERT:{0.insertion_side} {0.insertion_point}".format(self)

    @property
    def insert_details(self):
        return " ".join(map(str,(self.secondary_iv, self.secondary_iv.start,
                                 self.secondary_iv.end,
                                 self.secondary_iv.start_d,
                                 self.secondary_iv.end_d)))
    # def gap(self,other):
    #     return (other.iv.start - self.iv.end)-1
    def __iter__(self):
        return tuple(self.read_list)

    def next(self):
        return self.read_list.next()

    def draw(self):
        if self.insertion_side == LEFT:
            sorted_reads = sorted(self.read_list,reverse=True,key=lambda read: (read.iv.end,read.iv.start,))
        else:
            sorted_reads = sorted(self.read_list,key=lambda read: (read.iv.start,read.iv.end))
        return [read.draw(i) for i,read in enumerate(sorted_reads)]

class ClusterGroup:
    Header = ("# type",
              "ref_chrom",
              "left_start_b1", "left_junction_b1", "left_numreads",
              "gap_length", "insert_length",
              "insert_chrom", "insert_start", "insert_end", "insert_strand",
              "right_junction_b1","right_end_b1", "right_numreads")
    def __init__(self):
        pass

class ClusterDoublet(ClusterGroup):
    Type = "pair"
    def __init__(self,left_cluster, right_cluster):
        self.left = left_cluster
        self.right = right_cluster
        self.iv = left_cluster.iv.copy()
        self.iv.extend_to_include(right_cluster.iv)

        self.insert_iv = left_cluster.secondary_iv.copy()
        logger.debug("left_cluster:{0.secondary_iv}\tright_cluster:{1.secondary_iv}".format(left_cluster,right_cluster))
        self.insert_iv.extend_to_include(right_cluster.secondary_iv)
        if left_cluster.secondary_iv.start > right_cluster.secondary_iv.start:
            self.insert_iv.strand = "-"
        else:
            self.insert_iv.strand = "+"
        # self.gap = left_cluster.distance(right_cluster)

    def __str__(self):
        return "LEFT: {0.left}\nGAP: {0.gap_length} INSERT: {0.insert_length} {0.insert_iv} \nRIGHT: {0.right}\n\n".format(self)

    @property
    def str_with_secondary(self):
        return "LEFT: {0.left}, left_insert:{0.left.secondary_iv}\nGAP: {0.gap_length}\nRIGHT: {0.right}, right_insert: {0.right.secondary_iv}\n\n".format(self)

    @property
    def output(self):
        # ref_chrom, left_start_b1, left_junction_b1, gap_length, insert_length, insert_strand, insert_chrom, insert_start, insert_end, right_start, right_junction
        # print >>sys.stderr, "number of reads"; sys.exit(1)
        return (self.__class__.Type,
                self.left.iv.chrom,
                self.left.iv.start+1, self.left.iv.end,self.left.count,
                self.gap_length, self.insert_length,
                self.insert_iv.chrom, self.insert_iv.start+1, self.insert_iv.end,self.insert_iv.strand,
                self.right.iv.start+1, self.right.iv.end,self.right.count)

    @property
    def gap_length(self):
        return self.right.iv.start-self.left.iv.end

    @property
    def insert_length(self):
        return self.insert_iv.length
    
    def get_features(self,cluster_number):
        left_iv = self.left.iv.copy()
        left_iv.strand="+"
        left_feature = HTSeq.GenomicFeature("cluster_{0}_left".format(cluster_number), CLUSTER_GFF_TYPE, left_iv)
        left_feature.score = self.left.count
        
        right_iv = self.right.iv.copy()
        right_iv.strand="-"
        right_feature = HTSeq.GenomicFeature("cluster_{0}_right".format(cluster_number), CLUSTER_GFF_TYPE, right_iv)
        right_feature.score = self.right.count

        if self.gap == 0:
            gap_iv = HTSeq.GenomicInterval( left_iv.chrom, left_iv.end-1, right_iv.start+1, "." )
        else:
            # gap_iv = HTSeq.GenomicInterval( left_iv.chrom, left_iv.end+1, right_iv.start-1, "." )
            gap_iv = HTSeq.GenomicInterval( left_iv.chrom, left_iv.end, right_iv.start, "." )
        insert_name = "cluster_{0}_insert".format(cluster_number)
        insert_feature = HTSeq.GenomicFeature(insert_name, GAP_GFF_TYPE, gap_iv)
        # print dir(insert_feature)
        # insert_feature.__setattribute__("length",self.gap)
        insert_feature.attr = {"ID":insert_name, "length": self.gap}

        
        return (left_feature, insert_feature, right_feature)
        # attr:             The last (9th) column of a GFF file contains attributes, i.e. a list of name/value pairs. These are transformed into a dict, such that, e.g., gf.attr['gene_id'] gives the value of the attribute gene_id in the feature described by GenomicFeature object gf. The parser for the attribute field is reasonably flexible to deal with format variations (it was never clearly established whetehr name and value should be sperarated by a colon or an equal sign, and whether quotes need to be used) and also does a URL style decoding, as is often required.
    def draw(self):
        right_rect_list = self.right.draw()
        left_rect_list = self.left.draw()
        return right_rect_list+left_rect_list, max(len(right_rect_list),len(left_rect_list))
        
class ClusterSingleton(ClusterGroup):
    Type = "singleton"
    def __init__(self,singleton_cluster):
        self.singleton = singleton_cluster
        self.iv = singleton_cluster.iv

    def __str__(self):
        return "SINGLETON: {0.singleton}\n\n".format(self)

    @property
    def str_with_secondary(self):
        return "SINGLETON: {0.singleton}, insert:{0.singleton.secondary_iv}\n\n".format(self)

    def get_features(self,cluster_number):
        singleton_iv = self.singleton.iv.copy()
        insert_iv = None
        if self.singleton.insertion_side == RIGHT:
            singleton_iv.strand="+"
            # insert_iv = HTSeq.GenomicInterval(singleton_iv.chrom, singleton_iv.end+1, singleton_iv.end+1, "." )
            insert_iv = HTSeq.GenomicInterval(singleton_iv.chrom, singleton_iv.end, singleton_iv.end+1, "." )
        elif self.singleton.insertion_side == LEFT:
            singleton_iv.strand="-"
            # insert_iv =HTSeq.GenomicInterval(singleton_iv.chrom, singleton_iv.start-1, singleton_iv.start-1, "." )
            insert_iv =HTSeq.GenomicInterval(singleton_iv.chrom, singleton_iv.start-1, singleton_iv.start, "." )
        else:
            singleton_iv.strand="."

        singleton_feature = HTSeq.GenomicFeature("cluster_{0}_singleton".format(cluster_number), CLUSTER_GFF_TYPE, singleton_iv)
        singleton_feature.score = self.singleton.count
        feature_list = [singleton_feature] 

        if insert_iv != None:
            insert_feature = HTSeq.GenomicFeature("cluster_{0}_junction".format(cluster_number), SINGLETON_INSERT_GFF_TYPE, insert_iv)
            feature_list.append(insert_feature)
        
        return feature_list

    @property
    def output(self):
        if self.singleton.insertion_side == RIGHT:
            return (self.__class__.Type,
                    self.iv.chrom,
                    self.iv.start, self.iv.end,self.singleton.count,
                    None, None,
                    None,None,None,None,
                    None,None,None)
        elif self.singleton.insertion_side == LEFT:
            return (self.__class__.Type,
                    self.iv.chrom,
                    None,None,None,
                    None, None,
                    None,None,None,None,
                    self.iv.start, self.iv.end,self.singleton.count)
    def draw(self):
        rect_list = self.singleton.draw()
        return rect_list, len(rect_list)

if __name__ == "__main__":
   main()

