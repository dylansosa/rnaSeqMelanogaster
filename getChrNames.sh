for k in *assembly_report.txt
  do
    out=$(echo $k | sed 's/.txt/.chrnames/')
    grep -e '^[^#]' $k | awk '{ print $5, $7 }' > $out
done
