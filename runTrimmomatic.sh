for i in {57..73} ; do (f=$(echo SRR3305"$i"_1_*.fastq | sed s'/.fastq//') && r=$(echo SRR3305"$i"_2_*.fastq | sed s'/.fastq//') && trimmomatic PE  -threads 12 "$f".fastq "$r".fastq "$f".paired.fastq "$f".unpaired.fastq "$r".paired.fastq "$r".unpaired.fastq ILLUMINACLIP:adapters/illuminaUniversal.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36) ; done
