
454.quals.txt: ../../reads/454Junior-1.fastq ../../reads/454Junior-2.fastq
	cat ../../reads/454Junior-1.fastq ../../reads/454Junior-2.fastq | \
		extract_quals_from_fastq_by_position.py /dev/stdin "fastq" > 454.quals.txt

iontor.quals.txt: ../../reads/IonTorrent-1.fastq ../../reads/IonTorrent-2.fastq
	cat ../../reads/IonTorrent-1.fastq ../../reads/IonTorrent-2.fastq | \
		extract_quals_from_fastq_by_position.py /dev/stdin "fastq" > iontor.quals.txt

miseq.quals.txt: ../../reads/MiSeq-1.fastq
	cat ../../reads/MiSeq-1.fastq | \
		extract_quals_from_fastq_by_position.py /dev/stdin "fastq-illumina" > miseq.quals.txt

454.accuracy.txt: ../../alignments/454Junior-both.sorted.bam
	calculate_accuracy.py \
		../../alignments/454Junior-both.sorted.bam \
		../../assemblies/newbler/280_flxplus.fna \
		> 454.accuracy.txt

iontor.accuracy.txt: ../../alignments/IonTorrent-both.sorted.bam
	calculate_accuracy.py \
		../../alignments/IonTorrent-both.sorted.bam \
		../../assemblies/newbler/280_flxplus.fna \
		> iontor.accuracy.txt

miseq.accuracy.txt: calculate_accuracy.py \
                 ../../alignments/MiSeq-1_280.sorted.bam \
		 ../../assemblies/newbler/280_flxplus.fna \
		 > miseq.accuracy.txt


clean:
	rm 454.quals.txt 454_test.quals.txt iontor.quals.txt miseq.quals.txt


