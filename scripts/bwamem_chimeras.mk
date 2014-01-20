# make all  -l 16 -j # for parallel make
STAGE_DIR=staging
SEQ_DIR=seqs
# PZPNAT_SEQ=ref/pPZP-NATcc_plasmid_100813.fasta
PZPNAT_SEQ=$(COLLAB)/AlspaughLab/info/pPZP-NATcc.fasta 
H99_SEQ=$(STAGE_DIR)/cryptococcus_neoformans_grubii_h99_2_contigs.fasta.gz
MERGED_SEQ=$(SEQ_DIR)/h99_pzpnat.fa
INDEXBASE=$(basename $(MERGED_SEQ))
INDEX_FILES=$(addprefix $(INDEXBASE),.bwt .pac .ann .amb .sa)
FASTQ_DIR=/home/ske5/WHM_1/raw_fastqs
BWA_OUTDIR=bwa
NUMTHREADS=12
MINSEEDLEN=10
# BWA_OUTBAM=$(BWA_OUTDIR)/

.PRECIOUS: %.bam

dir_guard=@mkdir -p $(@D)

# all : test $(TOPHAT_BASE_DIR)/accepted_hits.bam $(STAGE_DIR)/ppzp-nat_reads.bam $(STAGE_DIR)/fusion_chrom_counts.txt
# all: $(BWA_OUTDIR)/Undetermined_S0.bam $(BWA_OUTDIR)/SE-WHM1_S1.bam
all: $(BWA_OUTDIR)/Undetermined_S0.unsrt.bam

test :
	echo $(INDEXBASE)
	echo $(INDEX_FILES)

# #===============================================================================
# # Download and merge reference genomes 
$(MERGED_SEQ) : $(H99_SEQ) $(PZPNAT_SEQ)
	$(dir_guard)
	zcat $(word 1,$^) | cat $(word 2,$^) - > $(STAGE_DIR)/$(@F)
	mv $(STAGE_DIR)/$(@F) $@

$(H99_SEQ) :
	$(dir_guard)
	# The following changes the chromosome names in the genome sequence file so they are compatible with the gtf file
	curl "http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EASupercontigsFasta&sp=SCNA2&sp=S.gz" | zcat | sed s/Supercontig_2/Chromosome_2/ | gzip -c > $@

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# Index reference genome
%.bwt %.pac %.ann %.amb %.sa : %.fa
	bwa index $< -p $*

#===============================================================================
# Run BWA
#-----------
# The following is used to make a comma separated list out of file list returned by 
# http://blog.jgc.org/2007/06/escaping-comma-and-space-in-gnu-make.html
, := ,
space :=
space +=

$(BWA_OUTDIR)/%.bam :  $(FASTQ_DIR)/%_L001_R1_001.fastq.gz $(FASTQ_DIR)/%_L001_R2_001.fastq.gz $(INDEX_FILES)
	echo $^
	$(dir_guard)
	$(eval TMPOUT := $(basename $@)_TMP)
	bwa mem -t $(NUMTHREADS) -k $(MINSEEDLEN) $(INDEXBASE) \
	$(word 1,$^) $(word 2,$^) | samtools view -Su - | samtools sort - $(TMPOUT)
	mv $(TMPOUT).bam $@

$(BWA_OUTDIR)/%.unsrt.bam :  $(FASTQ_DIR)/%_L001_R1_001.fastq.gz $(FASTQ_DIR)/%_L001_R2_001.fastq.gz $(INDEX_FILES)
	echo $^
	$(dir_guard)
	$(eval TMPOUT := $(basename $@)_TMP)
	bwa mem -t $(NUMTHREADS) -k $(MINSEEDLEN) $(INDEXBASE) \
	$(word 1,$^) $(word 2,$^) | samtools view -Sb - > $(TMPOUT)
	mv $(TMPOUT) $@
# #--------------------------------------------------------------------------------
# find supplementary alignments
# samtools view -f 0x800 Undetermined_S0_names.bam | grep PZP
# samtools view -f 0x40 Undetermined_S0_names.bam | grep M00830:103:000000000-A4A90:1:1103:25577:17255
# samtools view -f 0x80 Undetermined_S0_names.bam | grep M00830:103:000000000-A4A90:1:1103:25577:17255
# samtools view -f 0x80 Undetermined_S0_names.bam | grep M00830:103:000000000-A4A90:1:1103:25577:17255 | cut -f1-6,12-100
#===============================================================================
# Quick and dirty extract fusion reads
#-----------
$(STAGE_DIR)/ppzp-nat_reads.sam : $(TOPHAT_BASE_DIR)/accepted_hits.bam
	samtools view -H $< > $@
	samtools view $< | grep -i pPZP-NATcc >> $@

%.bam : %.sam
	samtools view -Su $< | samtools sort - $(basename $@)

#===============================================================================
# Quick and dirty fusion chrom counts
#-----------
$(STAGE_DIR)/fusion_chrom_counts.txt : %/accepted_hits.bam
	samtools view $< | grep "XF" | cut -f18 | cut -d" " -f2 | sort | uniq -c > $@

#===============================================================================
# Find fusion read junctions
#-----------
$(STAGE_DIR)/%_junctions.csv $(STAGE_DIR)/%_fusionreads.csv : %/accepted_hits.bam
	python2.7 $SCRIPTS/extract_reads_and_fusions.py $<  DUMMY --junction $*_junctions.csv --fusionreads $*_fusionreads.csv > /dev/null

# # Cleanup
# #-----------
# pristine:
# 	# rm -rf $(DIR_LIST)
# 	find . -type d -delete

# tidy:
# 	rm -rf $(STAGE_DIR)/*

# rm_th_results:
# 	rm -rf $(TOPHAT_BASE_DIR)/*
#================================================================================
#================================================================================
# Generate test data
# python2.7 /Users/josh/Documents/BioinfCollabs/scripts/make_random_sequence.py --output random_seqs.fa
# python2.7 /Users/josh/Documents/BioinfCollabs/scripts/make_synthetic_chimeric_reads.py random_seqs.fa --output random_seq_reads.fastq
# bwa index -p rand_seqs random_seqs.fa 
# bwa mem rand_seqs random_seq_reads.fastq | samtools view -Sb - > rand.bam
# python2.7 /Users/josh/Documents/BioinfCollabs/scripts/extract_bwa_chimeras.py rand.bam 
# python2.7 /Users/josh/Documents/BioinfCollabs/scripts/extract_bwa_chimeras.py rand.bam | egrep "prime|junc"
#--------------------------------------------------------------------------------
# python2.7 /Users/josh/Documents/BioinfCollabs/scripts/make_synthetic_chimeric_reads.py random_seqs.fa --output random_reads_3.fastq
# bwa mem rand_seqs random_reads_3.fastq | samtools view -Sb - > rand_3.bam
# python2.7 /Users/josh/Documents/BioinfCollabs/scripts/extract_bwa_chimeras.py rand_3.bam 
# python2.7 /Users/josh/Documents/BioinfCollabs/scripts/extract_bwa_chimeras.py rand_3.bam | egrep "prime|junc"
#--------------------------------------------------------------------------------