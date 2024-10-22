
FLOW_SRA_FILES=SRR388806.sra SRR388807.sra SRR389193.sra SRR389194.sra
FLOW_SFF_FILES=SRR388806.sff SRR388807.sff SRR389193.sff SRR389194.sff
FLOW_FASTQ_FILES=SRR388806.fastq SRR388807.fastq SRR389193.fastq SRR389194.fastq

ILLUMINA_SRA_FILES=SRR390202.sra
ILLUMINA_FASTQ_FILES=SRR390202_ACTTGA_1.fastq SRR390202_GATCAG_1.fastq SRR390202_TAGCTT_1.fastq SRR390202_CGATGT_1.fastq SRR390202_GGCTAC_1.fastq SRR390202_TTAGGC_1.fastq SRR390202_ACTTGA_2.fastq SRR390202_GATCAG_2.fastq SRR390202_TTAGGC_2.fastq SRR390202_ATCACG_2.fastq SRR390202_GGCTAC_2.fastq SRR390202_CGATGT_2.fastq SRR390202_TAGCTT_2.fastq
ILLUMINA_SHUFFLED_FILES=MiSeq-1_280.fastq MiSeq-1_282.fastq MiSeq-1_283.fastq MiSeq-1_518.fastq MiSeq-1_540.fastq MiSeq-1_541.fastq MiSeq-1_627.fastq
ILLUMINA_FILTERED_FILES=MiSeq-1_280_filtered.fastq

all: $(FLOW_FASTQ_FILES) $(ILLUMINA_FILTERED_FILES)

$(ILLUMINA_FILTERED_FILES): $(ILLUMINA_SHUFFLED_FILES)
	perl ../../deps/fastq_trim_trailing_Bs.pl < MiSeq-1_280.fastq > MiSeq-1_280_filtered.fastq

$(ILLUMINA_SHUFFLED_FILES): $(ILLUMINA_FASTQ_FILES)
	shuffle.py SRR390202_ATCACG_1.fastq SRR390202_ATCACG_2.fastq >MiSeq-1_280.fastq
	shuffle.py SRR390202_TTAGGC_1.fastq SRR390202_TTAGGC_2.fastq >MiSeq-1_282.fastq
	shuffle.py SRR390202_ACTTGA_1.fastq SRR390202_ACTTGA_2.fastq >MiSeq-1_283.fastq
	shuffle.py SRR390202_TAGCTT_1.fastq SRR390202_TAGCTT_2.fastq >MiSeq-1_518.fastq
	shuffle.py SRR390202_GATCAG_1.fastq SRR390202_GATCAG_2.fastq >MiSeq-1_540.fastq
	shuffle.py SRR390202_CGATGT_1.fastq SRR390202_CGATGT_2.fastq >MiSeq-1_541.fastq
	shuffle.py SRR390202_GGCTAC_1.fastq SRR390202_GGCTAC_2.fastq >MiSeq-1_627.fastq

$(FLOW_FASTQ_FILES): $(FLOW_SFF_FILES)

$(FLOW_SFF_FILES): $(FLOW_SRA_FILES)

$(ILLUMINA_FASTQ_FILES): $(ILLUMINA_SRA_FILES)
	fastq-dump -F --split-files -Q 64 -G --defline-seq "@\$sn#\$sg/\$ri" SRR390202.sra --defline-qual "+" $<

$(FLOW_SRA_FILES):
	get_from_sra.py $@

$(ILLUMINA_SRA_FILES):
	get_from_sra.py $@

%.sff: %.sra
	sff-dump $<

%.fastq: %.sff
	sff_extract_0_2_13 -c -Q $< -o $*

