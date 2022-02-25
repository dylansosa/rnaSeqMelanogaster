# for i in {57..73} ; do echo SRR3305"$i"_2_*.fastq ; done
#for i in {57..73} ; do echo SRR3305"$i"_1_*.fastq SRR3305"$i"_2_*.fastq paired_SRR3305"$i"_1_*.fastq unpaired_SRR3305"$i"_1_*.fastq paired_SRR3305"$i"_2_*.fastq unpaired_SRR3305"$i"_2_*.fastq ; donei

# for i in {58..58} ; do trimmomatic PE  -threads 12 SRR3305"$i"_1_*.fastq SRR3305"$i"_2_*.fastq paired_SRR3305"$i"_1_*.fastq unpaired_SRR3305"$i"_1_*.fastq paired_SRR3305"$i"_2_*.fastq unpaired_SRR3305"$i"_2_*.fastq ILLUMINACLIP:adapters/illuminaUniversal.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 ; done

# for i in {58..58} ; do name=$(echo SRR3305"$i"_1_*.fastq | cut -d '_' -f 2- | sed s'/.fastq//') ; done
# 1_Dpse_Female_Carcass_min_Reproductive_Tract_Replicate_2
# for i in {58..58} ; do (name=$(echo SRR3305"$i"_1_*.fastq | cut -d '_' -f 2- | sed s'/.fastq//') && trimmomatic PE  -threads 12 SRR3305"$i"_1_.fastq SRR3305"$i"_2_*.fastq paired_SRR3305"$i"_$name paired_SRR3305"$i"_2_*.fastq unpaired_SRR3305"$i"_2_*.fastq ILLUMINACLIP:adapters/illuminaUniversal.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36) ; done

#for i in {57..73} ; do (f=$(echo SRR3305"$i"_1_*.fastq | cut -d '_' -f 2- | sed s'/.fastq//') && r=$(echo SRR3305"$i"_2_*.fastq | cut -d '_' -f 2- | sed s'/.fastq//') && trimmomatic PE  -threads 12 "$f".fastq "$r".fastq "$f".paired.fastq "$f".unpaired.fastq "$r".paired.fastq "$r".unpaired.fastq ILLUMINACLIP:adapters/illuminaUniversal.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36) ; done

#for i in {57..73} ; do (f=$(echo SRR3305"$i"_1_*.fastq) && r=$(echo SRR3305"$i"_2_*.fastq) && trimmomatic PE  -threads 12 "$f".fastq "$r".fastq "$f".paired.fastq "$f".unpaired.fastq "$r".paired.fastq "$r".unpaired.fastq ILLUMINACLIP:adapters/illuminaUniversal.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36) ; done

for i in {57..73} ; do (f=$(echo SRR3305"$i"_1_*.fastq | sed s'/.fastq//') && r=$(echo SRR3305"$i"_2_*.fastq | sed s'/.fastq//') && trimmomatic PE  -threads 12 "$f".fastq "$r".fastq "$f".paired.fastq "$f".unpaired.fastq "$r".paired.fastq "$r".unpaired.fastq ILLUMINACLIP:adapters/illuminaUniversal.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36) ; done
