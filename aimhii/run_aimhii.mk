# make run_aimhii NUMTHREADS=12

##------------------------------------------------------------
## MISC
##------------------------------------------------------------
NUMTHREADS ?= 2
dir_guard=@mkdir -p $(@D)
RANDOM_SEED := 1
NUM_SUBSET ?= 200000
SRA_ACCESSION := SRR1964709

##------------------------------------------------------------
## GLOBAL DIRECTORIES
##------------------------------------------------------------
INFO_DIR := info
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
ADAPTER_FASTA :=$(INFO_DIR)/illumina_adapter1.fasta

READ1_FASTQ_SUBSET := $(FASTQ_DIR)/subset_$(SRA_ACCESSION)_R1.fastq.gz
READ2_FASTQ_SUBSET := $(FASTQ_DIR)/subset_$(SRA_ACCESSION)_R2.fastq.gz
SUBSET_OUTPUT := $(RESULTS_DIR)/$(notdir $(READ1_FASTQ_SUBSET:_R1.fastq.gz=_results.csv))

READ1_FASTQ := $(FASTQ_DIR)/full_$(SRA_ACCESSION)_R1.fastq.gz
READ2_FASTQ := $(FASTQ_DIR)/full_$(SRA_ACCESSION)_R2.fastq.gz
FULL_OUTPUT := $(RESULTS_DIR)/$(notdir $(READ1_FASTQ:_R1.fastq.gz=_results.csv))


## REFERENCE GENOME
##-----------------
PZPNAT_SEQ=$(INFO_DIR)/pPZP-NATcc.fasta
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
AIMHII ?= aimhii

##------------------------------------------------------------
## PHONY RULES
##------------------------------------------------------------
refgenome : $(H99_GENOME)

run_aimhii : $(FULL_OUTPUT)

run_subset : $(SUBSET_OUTPUT)

subset_data : $(READ1_FASTQ_SUBSET) $(READ2_FASTQ_SUBSET) 

full_data : $(READ1_FASTQ) $(READ2_FASTQ)

extract :
	extract_chimeras $(BWA_DIR)/merged.bam --insert $(PZPNAT_SEQ) --table $(RESULTS_DIR)/extract_results.csv 


#===============================================================================
# Download and merge reference genomes 
#===============================================================================
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
# $(RESULTS_DIR)/%_results.csv : $(BWA_DIR) $(H99_GENOME) $(FASTQ_DIR)/%_R1_001.fastq.gz $(FASTQ_DIR)/%_R2_001.fastq.gz 
$(RESULTS_DIR)/%_results.csv : $(BWA_DIR) $(H99_GENOME) $(FASTQ_DIR)/%_R1.fastq.gz $(FASTQ_DIR)/%_R2.fastq.gz 
	$(dir_guard)
	$(AIMHII) --threads $(NUMTHREADS) --outfile $@.tmp -t $(word 1,$^) $(word 2,$^) $(PZPNAT_SEQ) $(ADAPTER_FASTA) $(word 3,$^) $(word 4,$^)
	mv $@.tmp $@

#===============================================================================
# Generate sample subset
#===============================================================================
$(FASTQ_DIR)/subset_%_R1.fastq.gz $(FASTQ_DIR)/subset_%_R2.fastq.gz : 
	$(dir_guard)
	fastq-dump -X $(NUM_SUBSET) --split-files --gzip $* --outdir $(@D)
	mv $(@D)/$*_1.fastq.gz $(@D)/subset_$*_R1.fastq.gz
	mv $(@D)/$*_2.fastq.gz $(@D)/subset_$*_R2.fastq.gz

$(FASTQ_DIR)/full_%_R1.fastq.gz $(FASTQ_DIR)/full_%_R2.fastq.gz : 
	$(dir_guard)
	fastq-dump --split-files --gzip $* --outdir $(@D)
	mv $(@D)/$*_1.fastq.gz $(@D)/full_$*_R1.fastq.gz
	mv $(@D)/$*_2.fastq.gz $(@D)/full_$*_R2.fastq.gz

#===============================================================================
# CLEAN UP
#===============================================================================
clean : 
	rm -rf $(FASTQ_DIR) $(RESULTS_DIR) $(BWA_DIR)

pristine :
	rm -rf $(FASTQ_DIR) $(SEQ_DIR) $(RESULTS_DIR) $(BWA_DIR)

