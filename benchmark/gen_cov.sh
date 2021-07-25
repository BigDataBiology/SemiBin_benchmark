for i in $1/*_cov.txt
do
   echo $i
   stub=${i%_cov.txt}
   echo $stub
   awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > ${stub}_cov.csv
done
