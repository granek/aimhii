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
MAKE_INSERTION ?= make_synthetic_insertion
MAKE_RANDOM_SEQUENCE ?= make_random_sequence

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
	cat $@

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
# Synthetic Sequences
#===============================================================================
SCRIPTS_DIR=$(ALSPAUGH)/make_synthetic/scripts
SYN_DIR=synthetic
SYNTHETIC_FASTQ_DIR=$(SYN_DIR)
SYNTHETIC_FASTQS :=  syn_R1.fastq.gz syn_R2.fastq.gz
SYNTHETIC_FASTQ_FULLPATH := $(addprefix $(SYNTHETIC_FASTQ_DIR)/, $(SYNTHETIC_FASTQS))
SYNTHETIC_RESULTS := $(word 1, $(SYNTHETIC_FASTQ_FULLPATH:_R1.fastq.gz=_results.csv))
SYN_REFSEQ := $(SYN_DIR)/syn_ref.fa # a synthetic reference sequence
SYN_INSERTSEQ := $(SYN_DIR)/syn_insert.fa # a synthetic insert sequence
SYN_MUTSEQ := $(SYN_DIR)/syn_mut.fa # the  synthetic reference sequence with insertions of SYN_INSERTSEQ
SYNREF_AND_SYNINSERT_SEQ := $(SYN_DIR)/syn_concat.fa # concatenated synthetic reference and insert sequences (for mapping reference)
# SYNREF_AND_SYNMUT := $(SYN_DIR)/syn_ref_and_mut.fa # concatenated synthetic reference and SYN_MUTSEQ (for simulated read generation)
SYNREAD_LEN := 250
FASTQ_SUFFIX=.fastq.gz


syn : $(SYN_REFSEQ) $(SYN_INSERTSEQ) $(SYN_MUTSEQ) $(SYNTHETIC_FASTQ_FULLPATH) $(SYNTHETIC_RESULTS)
	echo $(SYNTHETIC_RESULTS)

$(SYN_DIR)/%_results.csv : $(BWA_DIR) $(SYN_REFSEQ) $(SYN_DIR)/%_R1.fastq.gz $(SYN_DIR)/%_R2.fastq.gz 
	$(dir_guard)
	$(AIMHII) --threads $(NUMTHREADS) --outfile $@.tmp -t $(word 1,$^) $(word 2,$^) $(SYN_INSERTSEQ) $(ADAPTER_FASTA) $(word 3,$^) $(word 4,$^)
	mv $@.tmp $@
	cat $@

%_R1$(FASTQ_SUFFIX) %_R2$(FASTQ_SUFFIX) : %_ref.fa %_mut.fa # $(SYN_REFSEQ) $(SYN_MUTSEQ)
	$(dir_guard)
	$(eval SEQ1_R1 := $(basename $(word 1,$^))_R1.fastq)
	$(eval SEQ1_R2 := $(basename $(word 1,$^))_R2.fastq)
	$(eval SEQ2_R1 := $(basename $(word 2,$^))_R1.fastq)
	$(eval SEQ2_R2 := $(basename $(word 2,$^))_R2.fastq)

	wgsim -S 1 -N 100000 -1 $(SYNREAD_LEN) -2 $(SYNREAD_LEN) $(word 1,$^) $(SEQ1_R1) $(SEQ1_R2)
	wgsim -S 2 -N 20000  -1 $(SYNREAD_LEN) -2 $(SYNREAD_LEN) $(word 2,$^) $(SEQ2_R1) $(SEQ2_R2)
	cat $(SEQ1_R1) $(SEQ2_R1) | gzip -c > $*_R1$(FASTQ_SUFFIX)
	cat $(SEQ1_R2) $(SEQ2_R2) | gzip -c > $*_R2$(FASTQ_SUFFIX)
	rm -f $(SEQ1_R1) $(SEQ1_R2) $(SEQ2_R1) $(SEQ2_R2)

$(SYNREF_AND_SYNINSERT_SEQ) : $(SYN_REFSEQ) $(SYN_INSERTSEQ)
	$(dir_guard)
	cat $(word 1,$^) $(word 2,$^) > $@.tmp
	mv $@.tmp $@

$(SYN_MUTSEQ) : $(SYN_REFSEQ) $(SYN_INSERTSEQ)
	$(dir_guard)
	$(MAKE_INSERTION) $(SYN_REFSEQ) $(SYN_INSERTSEQ) --position 90000 --delete 50 --output $@.tmp1
	$(MAKE_INSERTION) $(SYN_REFSEQ) $(SYN_INSERTSEQ) --position 30000 --delete 1 --output $@.tmp2
	$(MAKE_INSERTION) $(SYN_REFSEQ) $(SYN_INSERTSEQ) --position 10000 --delete 5 --invert --output $@.tmp3
	$(MAKE_INSERTION) $(SYN_REFSEQ)$(SYN_INSERTSEQ) --position 800 --delete 0 --output $@.tmp4
	cat $@.tmp1 $@.tmp2 $@.tmp3 $@.tmp4 > $@
	rm -f $@.tmp1 $@.tmp2 $@.tmp3 $@.tmp4

$(SYN_REFSEQ) :
	$(dir_guard)
	$(MAKE_RANDOM_SEQUENCE) --randseed 1 --prefix ref --lower --numseqs 1 --seqlen 100000 --output $@.tmp
	mv $@.tmp $@

$(SYN_INSERTSEQ) :
	$(dir_guard)
	$(MAKE_RANDOM_SEQUENCE) --randseed 2 --prefix insert --pad 10 --numseqs 1 --seqlen 2000 --output $@.tmp
	mv $@.tmp $@

cleansyn :
	rm -rf $(SYN_DIR)

#===============================================================================
# CLEAN UP
#===============================================================================
clean : 
	rm -rf $(FASTQ_DIR) $(RESULTS_DIR) $(BWA_DIR)

pristine :
	rm -rf $(FASTQ_DIR) $(SEQ_DIR) $(RESULTS_DIR) $(BWA_DIR)

