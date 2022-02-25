for i in *.fastq ; do species=$(echo $i | cut -f3 -d '_' | cut -c-3 ) && index=$(find ../transcriptomes/*$species* -type d -name *.idx) && SRR=$(echo $i | cut -f 1 -d '_') && salmon quant -i $index -l A -p 12 -r $i --validateMappings -o salmon/$SRR ; done
#--numBootstraps 100; done
