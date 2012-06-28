export PATH=$PATH:/mnt/fast/software/assembly/mira_3.4.0_prod_linux-gnu_x86_64_static/bin/

IonTorrent-1_assembly: parameters.txt
        ./sff_extract_0_2_10 -Q -s IonTorrent-1_in.iontor.fastq -x IonTorrent-1_traceinfo_in.iontor.xml ../../reads/IonTorrent-1.sff
        mira -project=IonTorrent-1 -params=parameters.txt >&log_assembly.txt

IonTorrent-2_assembly: parameters.txt
        ./sff_extract_0_2_10 -Q -s IonTorrent-2_in.iontor.fastq -x IonTorrent-2_traceinfo_in.iontor.xml ../../reads/IonTorrent-2.sff
        mira -project=IonTorrent-2 -params=parameters.txt >log_assembly.txt

IonTorrent-both_assembly: parameters.txt mira -project=IonTorrent-both -params=parameters.txt >log_assembly.txtt
        ./sff_extract_0_2_10 -Q -s IonTorrent-both_in.iontor.fastq -x IonTorrent-both_traceinfo_in.iontor.xml ../../reads/IonTorrent-1.sff ../../reads/IonTorrent-2.sff
        mira -project=IonTorrent-both -params=parameters.txt >log_assembly.txt

MiSeq-1_assembly: illumina_parameters.txt
        mira -project=MiSeq-1 -params=illumina_parameters.txt -s ../../reads/MiSeq-1_280.fastq >log_assembly.txt

454Junior-1_assembly: parameters.txt
        ./sff_extract_0_2_10 -Q -s 454Junior-1_in.454.fastq -x 454Junior-1_traceinfo_in.454.xml ../../reads/454Junior-1.sff
        mira -project=454Junior-1 -params=454_parameters.txt >log_assembly.txt
        rm -rf 454Junior-1_assembly/454Junior-1_d_tmp

454Junior-2_assembly: parameters.txt
        ./sff_extract_0_2_10 -Q -s 454Junior-2_in.454.fastq -x 454Junior-2_traceinfo_in.454.xml ../../reads/454Junior-2.sff
        mira -project=454Junior-2 -params=454_parameters.txt >log_assembly.txt
        rm -rf 454Junior-2_assembly/454Junior-2_d_tmp

454Junior-both_assembly: parameters.txt
        ./sff_extract_0_2_10 -Q -s 454Junior-both_in.454.fastq -x 454Junior-both_traceinfo_in.454.xml ../../reads/454Junior-1.sff ../../reads/454Junior-2.sff
        mira -project=454Junior-both -params=454_parameters.txt >log_assembly.txt
        rm -rf 454Junior-both_assembly/454Junior-both_d_tmp
