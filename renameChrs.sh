for i in *names ; do

echo "Renaming identifiers for: $i"

file="$(ls "$i" | cut -f 1 -d '.' | cut -f 2 -d '_')"
# echo $file 
# this is the part of the file name that is just numeric, using to match to assembly

assembly="$(find ../assemblies -type f -name *$file*.fna)"

cat $i | sed s'/^/sed -i s"\//' | sed s':\( [A-Z]\):\/\1:' | sed s'/\/ /\//' | sed s":$:\/\" $assembly:" > temp
# create file with each line as a sed command specific to the current assembly

chmod +x temp
# make executable the temp file that contains all of the sed commands

./temp
contigs="$(wc -l temp | cut -f 1 -d ' ')"
# run the substitution commands and save a variable to print the number of contigs renamed

echo "### $contigs Contigs Renamed ###" ; done
