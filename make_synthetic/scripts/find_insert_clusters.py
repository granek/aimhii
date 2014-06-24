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
    cluster_group_list = []
    ga = HTSeq.GenomicArray("auto", stranded=False,typecode="O")
    for cluster in cluster_list:
        ga[cluster.iv] = cluster
    cur_left = cur_gap = cur_right = None
    for cur_iv,cur_val in ga.steps():
        print >>sys.stderr, cur_iv, cur_val, cur_iv.length
        if cur_val == None:
            cur_gap = cur_iv
        elif cur_left == None:
            cur_left = cur_val
        else:
            cur_right = cur_val
            if cur_gap == None:
                new_group = ClusterDoublet(cur_left,cur_right)
                print >>sys.stderr, "Found a cluster doublet NO GAP:", 0, new_group.gap, new_group
                cur_left = cur_gap = cur_right = None
            elif cur_gap.length < max_gap_size:
                new_group = ClusterDoublet(cur_left,cur_right)
                print >>sys.stderr, "Found a cluster doublet:", cur_gap.length, new_group.gap, new_group
                cur_left = cur_gap = cur_right = None
            else:
                new_group = ClusterSingleton(cur_left)
                print >>sys.stderr, "Found a cluster singleton:", cur_gap.length, new_group
                cur_left = cur_right
                cur_right = cur_gap = None
            cluster_group_list.append(new_group)
    
    # cluster_list = sorted(cluster_list)
    # print "#"*80
    # for cluster in cluster_list:
    #     print "{0.count:<10} {0.range}".format(cluster)

    # # while (cluster_list):
    return cluster_group_list
        
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

def find_clusters(align_seq):
    cluster_list = []
    cur_cluster = None
    for aread in align_seq:
        if not aread.aligned:
            continue
        elif cur_cluster == None:
            #prime the pump
            cur_cluster = ReadCluster(aread)
            cluster_list.append(cur_cluster)
        elif cur_cluster.overlaps(aread):
            # aread overlaps cur_cluster, so it should be added
            cur_cluster.add_read(aread)
            # cur_cluster.add_read(aread)
            # cur_cluster_iv.extend_to_include(aread.iv)
        else:
            # aread doesn't overlap cur_cluster, so time to start a new cluster
            # print "Finished Cluster:", cur_cluster.range, cur_cluster.count
            cur_cluster = ReadCluster(aread)
            cluster_list.append(cur_cluster)
    # cluster_list.append(cur_cluster)
    # print cur_cluster.range, cur_cluster.count
    return cluster_list

class ReadCluster:
    def __init__(self,aread):
        self.read_list = [aread]
        # print aread
        self.iv = aread.iv.copy()
        self.iv.strand="."
        self._insert_point = aread.insert_point
        self._insert_side = aread.insert_side

    def add_read(self, aread):
        # print "self.iv", self.iv, type(self.iv)
        # print "aread", aread
        # print "aread.iv", aread.iv, type(aread.iv)
        self.iv.extend_to_include(aread.iv)
        self.read_list.append(aread)

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

    # def __lt__(self,other):
    #     if self.iv.chrom < 
    #     return self.iv < other.iv

    def distance(self,other):
        return max(self.iv.start-other.iv.end, other.iv.start-self.iv.end)

    def __str__(self):
        return "{0.iv.chrom}:{0.iv.start}-{0.iv.end} LENGTH:{0.iv.length} COUNT:{0.count} INSERT:{0.insertion_side} {0.insertion_point}".format(self)


class ClusterGroup:
    def __init__(self):
        pass

class ClusterDoublet(ClusterGroup):
    def __init__(self,left_cluster, right_cluster):
        self.left = left_cluster
        self.right = right_cluster
        self.gap = left_cluster.distance(right_cluster)

    def __str__(self):
        return "LEFT: {0}\nGAP: {1}\nRIGHT: {2}\n\n".format(self.left,self.gap,self.right)

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




        

        
class ClusterSingleton(ClusterGroup):
    def __init__(self,singleton_cluster):
        self.singleton = singleton_cluster

    def __str__(self):
        return "SINGLETON: {0}\n\n".format(self.singleton)

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


    
if __name__ == "__main__":
   main()

