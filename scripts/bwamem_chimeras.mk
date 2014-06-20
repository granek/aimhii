# /ssh:microbe:/home/josh/collabs/AlspaughLab/scripts/bwamem_chimeras.mk
# make all --load-average 8 --jobs --warn-undefined-variables # for parallel make
# make all  -l 16 -j # for parallel make
STAGE_DIR=staging
SEQ_DIR=seqs
BWA_OUTDIR=bwa
RESULTS_DIR=results

# PZPNAT_SEQ=ref/pPZP-NATcc_plasmid_100813.fasta
PZPNAT_SEQ=$(COLLAB)/AlspaughLab/info/pPZP-NATcc.fasta 
H99_SEQ=$(STAGE_DIR)/cryptococcus_neoformans_var._grubii_h99__cna3__3_supercontigs.fasta.gz
MERGED_SEQ=$(SEQ_DIR)/h99_pzpnat.fa
INDEXBASE=$(basename $(MERGED_SEQ))
INDEX_FILES=$(addprefix $(INDEXBASE),.bwt .pac .ann .amb .sa)
RAW_FASTQ_DIR=/home/ske5/WHM_1/raw_fastqs
FASTQ_DIR= $(RAW_FASTQ_DIR)
# TRIM_FASTQ_DIR=trim_fastqs
##----------------------------------------------------------------------
FINAL_FASTQ_DIR=final_fastqs
ADAPTER_FASTA=info/illumina_adapter1.fasta
#--------------------------------------------------
# RAW_FASTQS_FULLPATH := 
#--------------------------------------------------
# /home/ske5/WHM_1/raw_fastqs/SE-WHM1_S1_L001_R1_001.fastq.gz
# /home/ske5/WHM_1/raw_fastqs/SE-WHM1_S1_L001_R2_001.fastq.gz
# /home/ske5/WHM_1/raw_fastqs/Undetermined_S0_L001_R1_001.fastq.gz
# /home/ske5/WHM_1/raw_fastqs/Undetermined_S0_L001_R2_001.fastq.gz


ALL_FASTQS := $(notdir $(wildcard $(RAW_FASTQ_DIR)/*.fastq.gz))
UNDETERMINED_FASTQS := $(notdir $(wildcard $(RAW_FASTQ_DIR)/Undetermined*.fastq.gz))
FASTQS := $(UNDETERMINED_FASTQS)
FASTQ_SUFFIX=_001.fastq.gz

TRIM_FASTQS := $(patsubst %$(FASTQ_SUFFIX),$(FINAL_FASTQ_DIR)/%.trim.fastq.gz,$(FASTQS))
R1_TRIM_FASTQS := $(patsubst %_R1.trim.fastq.gz,%,$(filter %_R1.trim.fastq.gz, $(TRIM_FASTQS)))
JOIN_FASTQS := $(addsuffix .un1.fastq.gz , $(R1_TRIM_FASTQS)) $(addsuffix .un2.fastq.gz ,$(R1_TRIM_FASTQS)) $(addsuffix .join.fastq.gz ,$(R1_TRIM_FASTQS))
UNSRT_BAMS := $(addprefix $(BWA_OUTDIR)/, $(addsuffix .join.unsrtbam,$(notdir $(R1_TRIM_FASTQS))) $(addsuffix .pair.unsrtbam,$(notdir $(R1_TRIM_FASTQS))))
SORTED_BAMS := $(addprefix $(BWA_OUTDIR)/, $(addsuffix .join.bam,$(notdir $(R1_TRIM_FASTQS))) $(addsuffix .pair.bam,$(notdir $(R1_TRIM_FASTQS))))
MERGED_BAMS := $(addprefix $(BWA_OUTDIR)/, $(addsuffix .merge.bam,$(notdir $(R1_TRIM_FASTQS))))
MERGE_JUNCS = $(addprefix $(RESULTS_DIR)/,  $(subst .bam,.junc,$(notdir $(MERGED_BAMS))))


# # R1_TRIM_FASTQS := $(patsubst %R1.trim.fastq.gz,%,$(TRIM_FASTQS))
# JOIN_FASTQS := $(addsuffix .un1.fastq.gz , $(R1_TRIM_FASTQS)) $(addsuffix .un2.fastq.gz ,$(R1_TRIM_FASTQS)) $(addsuffix .join.fastq.gz ,$(R1_TRIM_FASTQS))
# UNSRT_BAMS := $(addprefix  $(BWA_OUTDIR)/, $(addsuffix .join.bam,$(notdir $(R1_TRIM_FASTQS))) $(addsuffix .pair.bam,$(notdir $(R1_TRIM_FASTQS))))
# PZP_BAMS := $(patsubst %.bam,%.PZP.bam,$(UNSRT_BAMS))
# SA_BAMS := $(patsubst %.bam,%.SA.bam,$(UNSRT_BAMS))
# PZP_SA_BAMS := $(patsubst %.bam,%.PZP_SA.bam,$(UNSRT_BAMS))
# PZP_SA_SORT_BAMS := $(patsubst %.bam,%.sort.bam,$(PZP_SA_BAMS))
# PZP_SA_SORT_BAIS := $(patsubst %.bam,%.bam.bai,$(PZP_SA_SORT_BAMS))
# PZP_MERGE = $(BWA_OUTDIR)/SE-WHM1.PZP.merge.bam
# PZP_SA_MERGE = $(BWA_OUTDIR)/SE-WHM1.PZP_SA.merge.bam
# PZP_MERGE_SORT = $(BWA_OUTDIR)/SE-WHM1.PZP.merge.sort.bam
# PZP_SA_MERGE_SORT = $(BWA_OUTDIR)/SE-WHM1.PZP_SA.merge.sort.bam
# SA_MERGE_SORT = $(BWA_OUTDIR)/SE-WHM1.SA.merge.sort.bam
# PZP_SA_MERGE_JUNCS = $(RESULTS_DIR)/SE-WHM1.PZP_SA.merge.junc
# MERGE_JUNCS = $(RESULTS_DIR)/SE-WHM1.merge.junc
# JOIN_FASTQS := $(addprefix $(TRIM_FASTQS),.un1.fastq.gz .un2.fastq.gz .join.fastq.gz)


# # TRIM_FASTQS := $(addprefix $(FINAL_FASTQ_DIR)/,$(FASTQS))
# #--------------------------------------------------
# TOPHAT_BASE_DIR=thout

# FINAL_BAMS := $(patsubst %$(FASTQ_SUFFIX),$(TOPHAT_BASE_DIR)/%/accepted_hits.bam,$(FASTQS))
# FINAL_BAIS := $(addsuffix .bai,$(FINAL_BAMS))
# ##----------------------------------------------------------------------



NUMTHREADS=12
MINSEEDLEN=10
# BWA_OUTBAM=$(BWA_OUTDIR)/

.PRECIOUS: %.bam %.un1.fastq.gz %.un2.fastq.gz %.join.fastq.gz %.bwt %.pac %.ann %.amb %.sa

dir_guard=@mkdir -p $(@D)

# all : test $(TOPHAT_BASE_DIR)/accepted_hits.bam $(STAGE_DIR)/ppzp-nat_reads.bam $(STAGE_DIR)/fusion_chrom_counts.txt
# all: $(BWA_OUTDIR)/Undetermined_S0.bam $(BWA_OUTDIR)/SE-WHM1_S1.bam
all: todo $(TRIM_FASTQS) $(JOIN_FASTQS) $(UNSRT_BAMS) $(SORTED_BAMS) $(MERGED_BAMS) $(MERGE_JUNCS) # $(PZP_MERGE_SORT) $(PZP_SA_MERGE_SORT) $(PZP_SA_MERGE_JUNCS) $(MERGE_JUNCS) $(SA_BAMS) $(SA_MERGE_SORT) # $(BWA_OUTDIR)/Undetermined_S0.unsrt.bam 

todo :
	@echo "Get sequences from NCBI instead of Broad for long term stability"
	@echo "Test that sequence was actually downloaded.  Maybe try MD5 http://www.kolpackov.net/pipermail/notes/2004-September/000011.html"
test :
	@echo "-----INDEXBASE-----"
	@echo $(INDEXBASE)
	@echo "-----INDEX_FILES-----"
	@echo $(INDEX_FILES)
	@echo "-----JOIN_FASTQS-----"
	@echo $(JOIN_FASTQS)
	@echo "-----R1_TRIM_FASTQS-----"
	@echo $(R1_TRIM_FASTQS)
	@echo "-----UNSRT_BAMS-----"
	@echo $(UNSRT_BAMS)


# #===============================================================================
# # Download and merge reference genomes 
$(MERGED_SEQ) : $(H99_SEQ) $(PZPNAT_SEQ)
	$(dir_guard)
	zcat $(word 1,$^) | cat $(word 2,$^) - > $(STAGE_DIR)/$(@F)
	mv $(STAGE_DIR)/$(@F) $@



H99_SEQ_URL := "http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EASupercontigsFasta&sp=SCNA3&sp=S.gz"
GTF_URL := "http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EATranscriptsGtf&sp=SCNA3&sp=S.gz"
$(H99_SEQ) :
	$(dir_guard)
	# The following changes the chromosome names in the genome sequence file so they are compatible with the gtf file
	## curl "http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EASupercontigsFasta&sp=SCNA2&sp=S.gz" | zcat | sed s/Supercontig_2/Chromosome_2/ | gzip -c > $@
	curl -o $@.tmp $(H99_SEQ_URL) 
	mv $@.tmp $@

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



$(BWA_OUTDIR)/%.join.unsrtbam :  $(FINAL_FASTQ_DIR)/%.join.fastq.gz $(INDEX_FILES)
	# echo $^
	$(dir_guard)
	$(eval TMPOUT := $(basename $@)_TMP)
	bwa mem -t $(NUMTHREADS) -k $(MINSEEDLEN) $(INDEXBASE) \
	$(word 1,$^) | samtools view -Sb - > $(TMPOUT)
	mv $(TMPOUT) $@

$(BWA_OUTDIR)/%.pair.unsrtbam :  $(FINAL_FASTQ_DIR)/%.un1.fastq.gz $(FINAL_FASTQ_DIR)/%.un2.fastq.gz $(INDEX_FILES)
	# echo $^
	$(dir_guard)
	$(eval TMPOUT := $(basename $@)_TMP)
	bwa mem -t $(NUMTHREADS) -k $(MINSEEDLEN) $(INDEXBASE) \
	$(word 1,$^) $(word 2,$^) | samtools view -Sb - > $(TMPOUT)
	mv $(TMPOUT) $@

%.bam : %.unsrtbam
	$(dir_guard)
	$(eval TMPOUT := $(basename $@)_TMP)
	samtools sort $<  $(TMPOUT)
	mv $(TMPOUT).bam $@

# %.SA.sam : %.bam
# 	$(dir_guard)
# 	samtools view -H $< > $@.tmp
# 	samtools view $< | grep "SA:Z" >> $@.tmp
# 	mv $@.tmp $@

# %.PZP_SA.sam : %.PZP.bam
# 	$(dir_guard)
# 	samtools view -H $< > $@.tmp
# 	samtools view $< | grep "SA:Z" >> $@.tmp
# 	mv $@.tmp $@

# %.PZP.sam : %.bam
# 	$(dir_guard)
# 	samtools view -H $< > $@.tmp
# 	samtools view $< | grep "PZP" >> $@.tmp
# 	mv $@.tmp $@

# %.bam : %.sam
# 	$(dir_guard)
# 	samtools view -Sb > $@.tmp
# 	mv $@.tmp $@



# %.sort.bam : %.bam
# 	$(dir_guard)
# 	$(eval TMPOUT := $(basename $@)_TMP)
# 	samtools sort $<  $(TMPOUT)
# 	mv $(TMPOUT).bam $@

%.bam.bai : %.bam
	samtools index $<

%.merge.bam : %.pair.bam %.join.bam
	samtools merge $@ $^



# $(PZP_SA_MERGE) : $(PZP_SA_BAMS)
# 	samtools merge $@ $^

# $(PZP_MERGE) : $(PZP_BAMS)
# 	samtools merge $@ $^

# $(SA_MERGE) : $(SA_BAMS)
# 	samtools merge $@ $^
#----------------------------------------

# $(BWA_OUTDIR)/%.bam :  $(FASTQ_DIR)/%_L001_R1_001.fastq.gz $(FASTQ_DIR)/%_L001_R2_001.fastq.gz $(INDEX_FILES)
# 	echo $^
# 	$(dir_guard)
# 	$(eval TMPOUT := $(basename $@)_TMP)
# 	bwa mem -t $(NUMTHREADS) -k $(MINSEEDLEN) $(INDEXBASE) \
# 	$(word 1,$^) $(word 2,$^) | samtools view -Su - | samtools sort - $(TMPOUT)
# 	mv $(TMPOUT).bam $@

# $(BWA_OUTDIR)/%.unsrt.bam :  $(FASTQ_DIR)/%_L001_R1_001.fastq.gz $(FASTQ_DIR)/%_L001_R2_001.fastq.gz $(INDEX_FILES)
# 	echo $^
# 	$(dir_guard)
# 	$(eval TMPOUT := $(basename $@)_TMP)
# 	bwa mem -t $(NUMTHREADS) -k $(MINSEEDLEN) $(INDEXBASE) \
# 	$(word 1,$^) $(word 2,$^) | samtools view -Sb - > $(TMPOUT)
# 	mv $(TMPOUT) $@
# #--------------------------------------------------------------------------------
# find supplementary alignments
# samtools view -f 0x800 Undetermined_S0_names.bam | grep PZP
# samtools view -f 0x40 Undetermined_S0_names.bam | grep M00830:103:000000000-A4A90:1:1103:25577:17255
# samtools view -f 0x80 Undetermined_S0_names.bam | grep M00830:103:000000000-A4A90:1:1103:25577:17255
# samtools view -f 0x80 Undetermined_S0_names.bam | grep M00830:103:000000000-A4A90:1:1103:25577:17255 | cut -f1-6,12-100
#===============================================================================
# # Quick and dirty extract fusion reads
# #-----------
# $(STAGE_DIR)/ppzp-nat_reads.sam : $(TOPHAT_BASE_DIR)/accepted_hits.bam
# 	samtools view -H $< > $@
# 	samtools view $< | grep -i pPZP-NATcc >> $@

# %.bam : %.sam
# 	samtools view -Su $< | samtools sort - $(basename $@)

# #===============================================================================
# # Quick and dirty fusion chrom counts
# #-----------
# $(STAGE_DIR)/fusion_chrom_counts.txt : %/accepted_hits.bam
# 	samtools view $< | grep "XF" | cut -f18 | cut -d" " -f2 | sort | uniq -c > $@

# #===============================================================================
# # Find fusion read junctions
# #-----------
# $(STAGE_DIR)/%_junctions.csv $(STAGE_DIR)/%_fusionreads.csv : %/accepted_hits.bam
# 	python2.7 $(SCRIPTS)/extract_reads_and_fusions.py $<  DUMMY --junction $*_junctions.csv --fusionreads $*_fusionreads.csv > /dev/null

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
#--------------------------------------------------------------------------------
# Trim adapters
$(FINAL_FASTQ_DIR)/%_R1.trim.fastq.gz $(FINAL_FASTQ_DIR)/%_R2.trim.fastq.gz : $(ADAPTER_FASTA) $(RAW_FASTQ_DIR)/%_R1$(FASTQ_SUFFIX) $(RAW_FASTQ_DIR)/%_R2$(FASTQ_SUFFIX)
	$(dir_guard)
	fastq-mcf $< $(word 2,$^) $(word 3,$^) -o $(FINAL_FASTQ_DIR)/$*_R1.tmp.gz -o $(FINAL_FASTQ_DIR)/$*_R2.tmp.gz
	mv $(FINAL_FASTQ_DIR)/$*_R1.tmp.gz $(FINAL_FASTQ_DIR)/$*_R1.trim.fastq.gz
	mv $(FINAL_FASTQ_DIR)/$*_R2.tmp.gz $(FINAL_FASTQ_DIR)/$*_R2.trim.fastq.gz
#--------------------------------------------------------------------------------
# End join pairs
$(FINAL_FASTQ_DIR)/%.un1.fastq.gz $(FINAL_FASTQ_DIR)/%.un2.fastq.gz $(FINAL_FASTQ_DIR)/%.join.fastq.gz : $(FINAL_FASTQ_DIR)/%_R1.trim.fastq.gz $(FINAL_FASTQ_DIR)/%_R2.trim.fastq.gz
	$(dir_guard)
	fastq-join $(word 1,$^) $(word 2,$^) -o $(FINAL_FASTQ_DIR)/$*.%.tmp.gz
	mv $(FINAL_FASTQ_DIR)/$*.un1.tmp.gz $(FINAL_FASTQ_DIR)/$*.un1.fastq.gz
	mv $(FINAL_FASTQ_DIR)/$*.un2.tmp.gz $(FINAL_FASTQ_DIR)/$*.un2.fastq.gz
	mv $(FINAL_FASTQ_DIR)/$*.join.tmp.gz $(FINAL_FASTQ_DIR)/$*.join.fastq.gz
#--------------------------------------------------------------------------------
## python2.7 /home/josh/collabs/scripts/extract_bwa_chimeras.py Undetermined_S0_L001.pair.bam --junction Undetermined_S0_L001.pair.junc --fusionreads Undetermined_S0_L001.pair.fusread
#-----------
$(RESULTS_DIR)/%.junc $(STAGE_DIR)/%.fusread : $(BWA_OUTDIR)/%.bam
	$(dir_guard)
	# python2.7 $SCRIPTS/extract_reads_and_fusions.py $<  DUMMY --junction $*_junctions.csv --fusionreads $*_fusionreads.csv > /dev/null
	python2.7 $(SCRIPTS)/extract_bwa_chimeras.py $< --junction $(dir $@)$*.junc.tmp --fusionreads $(dir $@)$*.fusread.tmp
	mv $(dir $@)$*.junc.tmp $(dir $@)$*.junc
	mv $(dir $@)$*.fusread.tmp $(dir $@)$*.fusread


reset : 
	rm -rf $(FINAL_FASTQ_DIR) $(STAGE_DIR) $(SEQ_DIR) $(BWA_OUTDIR) $(RESULTS_DIR)

cleancurrent :
	rm -rf $(STAGE_DIR) $(SEQ_DIR) $(BWA_OUTDIR) $(RESULTS_DIR)
