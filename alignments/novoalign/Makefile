
280_flxplus.ndx: ../assemblies/newbler/280_flxplus.fna
	novoindex -k 14 -s 1 280_flxplus.ndx ../assemblies/newbler/280_flxplus.fna

454Junior-both.sam: 280_flxplus.ndx ../../analysis/fastqc/454Junior-both.fastq
	novoalign -d 280_flxplus.ndx -f ../../analysis/fastqc/454Junior-both.fastq -g 20 -x 5 -n 300 \
		-o SAM -K 454Junior-both.qcal.csv --hpstats 454Junior-both.indel.csv \
		> 454Junior-1.sam

IonTorrent-both.sam: 280_flxplus.ndx ../../analysis/fastqc/IonTorrent-both.fastq
	novoalign -d 280_flxplus.ndx -f ../../analysis/fastqc/IonTorrent-both.fastq -g 20 -x 5 -n 300 \
		-o SAM -K IonTorrent-both.qcal.csv --hpstats IonTorrent-both.indel.csv \
		> IonTorrent-both.sam

all: 454Junior-both.sam IonTorrent-both.sam

