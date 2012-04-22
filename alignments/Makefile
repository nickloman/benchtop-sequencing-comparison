
454Junior-both.sorted.bam: 454Junior-1.sorted.bam 454Junior-2.sorted.bam
	samtools merge 454Junior-both.sorted.bam 454Junior-1.sorted.bam 454Junior-2.sorted.bam
	samtools index 454Junior-both.sorted.bam

454Junior-1.sorted.bam: ../assemblies/newbler/280_flxplus.fna.fai
	gobwa 454Junior-1 ../reads/454Junior-1.fastq ../assemblies/newbler/280_flxplus.fna

454Junior-2.sorted.bam: ../assemblies/newbler/280_flxplus.fna.fai
	gobwa 454Junior-2 ../reads/454Junior-1.fastq ../assemblies/newbler/280_flxplus.fna

IonTorrent-both.sorted.bam: IonTorrent-1.sorted.bam IonTorrent-2.sorted.bam
	samtools merge IonTorrent-both.sorted.bam IonTorrent-1.sorted.bam IonTorrent-2.sorted.bam
	samtools index IonTorrent-both.sorted.bam

IonTorrent-1.sorted.bam: ../assemblies/newbler/280_flxplus.fna.fai
	gobwa IonTorrent-1 ../reads/IonTorrent-1.fastq ../assemblies/newbler/280_flxplus.fna

IonTorrent-2.sorted.bam: ../assemblies/newbler/280_flxplus.fna.fai
	gobwa IonTorrent-2 ../reads/IonTorrent-2.fastq ../assemblies/newbler/280_flxplus.fna

MiSeq-1_280.sorted.bam: ../assemblies/newbler/280_flxplus.fna.fai
	gobwa MiSeq-1_280 ../reads/MiSeq-1_280.fastq ../assemblies/newbler/280_flxplus.fna

../assemblies/newbler/280_flxplus.fna.fai: ../assemblies/newbler/280_flxplus.fna
	bwa index -a is ../assemblies/newbler/280_flxplus.fna
	samtools faidx ../assemblies/newbler/280_flxplus.fna

clean:
	rm 454Junior-both.sorted.bam 454Junior-1.sorted.bam 454Junior-2.sorted.bam

all: 454Junior-both.sorted.bam IonTorrent-both.sorted.bam MiSeq-1_280.sorted.bam
