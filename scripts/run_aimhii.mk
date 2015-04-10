# /ssh:gems:/home/josh/collabs/AlspaughLab/scripts/run_aimhii.mk
# make run_aimhii NUMTHREADS=12



##------------------------------------------------------------
## MISC
##------------------------------------------------------------
NUMTHREADS ?= 12
dir_guard=@mkdir -p $(@D)
RANDOM_SEED := 1
NUM_SUBSET := 400000

##------------------------------------------------------------
## GLOBAL DIRECTORIES
##------------------------------------------------------------
BASE_DIR := $(COLLAB)/AlspaughLab
AIMHII_DIR := $(BASE_DIR)/make_synthetic/scripts
ORIGINAL_FASTQ_DIR := /nfs/gems_sata/alspaugh/raw_fastqs
FASTQ_DIR := raw_fastqs
##------------------------------------------------------------
## LOCAL DIRECTORIES
##------------------------------------------------------------
TEMP_DIR := tempdir
SEQ_DIR := genome
RESULTS_DIR := results
BWA_DIR := bwa


##------------------------------------------------------------
## FILES
##------------------------------------------------------------
AIMHII_OUTPUT := $(RESULTS_DIR)/full_aimhii_results.csv
ADAPTER_FASTA :=$(BASE_DIR)/info/illumina_adapter1.fasta
READ1_FASTQ := $(FASTQ_DIR)/SE-WHM1_and_Undetermined_R1_001.fastq.gz
READ2_FASTQ := $(FASTQ_DIR)/SE-WHM1_and_Undetermined_R2_001.fastq.gz
FULL_OUTPUT := $(RESULTS_DIR)/$(notdir $(READ1_FASTQ:_R1_001.fastq.gz=_results.csv))


READ1_FASTQ_SUBSET := $(FASTQ_DIR)/sample_R1_001.fastq.gz
READ2_FASTQ_SUBSET := $(FASTQ_DIR)/sample_R2_001.fastq.gz
SUBSET_OUTPUT := $(RESULTS_DIR)/$(notdir $(READ1_FASTQ_SUBSET:_R1_001.fastq.gz=_results.csv))


## REFERENCE GENOME
##-----------------
PZPNAT_SEQ=$(BASE_DIR)/info/pPZP-NATcc.fasta
H99_SEQNAME := GCF_000149245.1_CNA3_genomic.fna.gz
H99_GENOME := $(SEQ_DIR)/$(basename $(H99_SEQNAME))
H99_SEQ=$(TEMP_DIR)/$(H99_SEQNAME)
H99_AND_PZPNAT_SEQ := $(SEQ_DIR)/h99_pzpnat.fa


# http://www.ncbi.nlm.nih.gov/genome?LinkName=pubmed_genome&from_uid=24743168
## rsync --partial --progress -av ftp.ncbi.nlm.nih.gov::genomes/all/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.fna.gz 
# H99_SEQ_URL := "http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EASupercontigsFasta&sp=SCNA3&sp=S.gz"
# GTF_URL := "http://www.broadinstitute.org/annotation/genome/cryptococcus_neoformans/download/?sp=EATranscriptsGtf&sp=SCNA3&sp=S.gz"
H99_SEQ_URL := "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000149245.1_CNA3/$(H99_SEQNAME)"
# GFF_URL := "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF_000149245.1_CNA3/GCF_000149245.1_CNA3_genomic.gff.gz"


##------------------------------------------------------------
## BINARIES
##------------------------------------------------------------
AIMHII := $(AIMHII_DIR)/aimhii.py 


##------------------------------------------------------------
## PHONY RULES
##------------------------------------------------------------
refgenome : $(H99_GENOME)

run_aimhii : $(FULL_OUTPUT)

test : $(SUBSET_OUTPUT)

subset : $(READ1_FASTQ_SUBSET) $(READ2_FASTQ_SUBSET) 

check :
	echo $(ORIGINAL_FASTQ_DIR)

extract :
	$(AIMHII_DIR)/extract_bwa_chimeras.py $(BWA_DIR)/merged.bam --insert $(PZPNAT_SEQ) --table $(RESULTS_DIR)/extract_results.csv 


#===============================================================================
# Download and merge reference genomes 
#===============================================================================
# $(H99_AND_PZPNAT_SEQ) : $(H99_SEQ) $(PZPNAT_SEQ)
# 	$(dir_guard)
# 	zcat $(word 1,$^) | cat $(word 2,$^) - > $(TEMP_DIR)/$(@F)
# 	mv $(TEMP_DIR)/$(@F) $@

# $(H99_SEQ) :
# 	$(dir_guard)
# 	curl --fail -o $@.tmp $(H99_SEQ_URL) 
# 	mv $@.tmp $@

$(H99_GENOME) :
	$(dir_guard)
	curl --fail -o $@.tmp.gz $(H99_SEQ_URL)
	gunzip $@.tmp.gz
	mv $@.tmp $@

#===============================================================================
# Setup
#===============================================================================
$(BWA_DIR) :
	mkdir $@

#===============================================================================
# Run It!
#===============================================================================
# $(AIMHII_OUTPUT) : $(H99_GENOME) $(BWA_DIR)
# $(AIMHII_OUTPUT) : $(H99_GENOME) $(BWA_DIR)
$(RESULTS_DIR)/%_results.csv : $(BWA_DIR) $(H99_GENOME) $(FASTQ_DIR)/%_R1_001.fastq.gz $(FASTQ_DIR)/%_R2_001.fastq.gz 
	$(dir_guard)
	$(AIMHII) --threads $(NUMTHREADS) --outfile $@.tmp -t $(word 1,$^) $(word 2,$^) $(PZPNAT_SEQ) $(ADAPTER_FASTA) $(word 3,$^) $(word 4,$^)
	mv $@.tmp $@

#===============================================================================
# Generate sample subset
#===============================================================================
$(FASTQ_DIR)/sample_R%_001.fastq.gz : $(FASTQ_DIR)/SE-WHM1_and_Undetermined_R%_001.fastq.gz
	$(dir_guard)
	# seqtk sample -s $(RANDOM_SEED) $< $(NUM_SUBSET) > $(basename $@).tmp
	# gzip  $(basename $@).tmp
	# mv $(basename $@).tmp.gz $@
	seqtk sample -s $(RANDOM_SEED) $< $(NUM_SUBSET) | gzip > $@.tmp.gz
	mv $@.tmp.gz $@


$(FASTQ_DIR)/SE-WHM1_and_Undetermined_R%_001.fastq.gz : $(ORIGINAL_FASTQ_DIR)/SE-WHM1_and_Undetermined_R%_001.fastq.gz
	$(dir_guard)
	ln -s $< $@

#===============================================================================
# CLEAN UP
#===============================================================================
clean : 
	rm -rf $(FASTQ_DIR) $(RESULTS_DIR) $(BWA_DIR)

pristine :
	rm -rf $(FASTQ_DIR) $(SEQ_DIR) $(RESULTS_DIR) $(BWA_DIR)

