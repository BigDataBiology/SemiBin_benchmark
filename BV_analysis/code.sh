# by Chengkaizhu
# 07282021

########################################################
# fastani
########################################################
# input fna from bins
fastANI --ql $file --rl $file -o $out --fragLen 1000 -t 80 --matrix

########################################################
# prokka
########################################################
# input fa 
for i in fasta/*.fa;do
    file=${i##*/}
    base=${file%.fna}
    #echo $base
    prokka $i --outdir prokka/$base --prefix $base --metagenome --cpus 20 --force --kingdom Bacteria
done
########################################################
# roary
########################################################
# input prokaa gff results 
	
# crc
roary *.gff --mafft -p 64 -i 95 -cd 100 -f $roary/roary_crc_out -e
query_pan_genome -a union $gff -o $roary/roary_crc_out/union.txt
# GMGC
roary *.gff --mafft -p 64 -i 95 -cd 99 -f $roary/roary_GMGC_out -e -g 150000
query_pan_genome -a union $gff -o $roary/roary_GMGC_crc/union.txt

########################################################
# scoary
########################################################
# input roary gene_presence_absence.csv result

scoary -g $roary/roary_crc_out/gene_presence_absence.csv -t human1dog0.csv -c BH -o $scoary/significant_genes

########################################################
# Iqtree
########################################################
# input roary core_gene result

iqtree -s $roary/roary_crc_out/core_gene_alignment.aln -m MF -nt AUTO

Best_fit=$(cat $roary/roary_GMGC_out/core_gene_alignment.aln.log|grep "Best-fit model"|awk '{print $3}')

echo "Best fit model is ${Best_fit}"

# GTR+F+I+G4
iqtree -s $roary/roary_crc_out/core_gene_alignment.aln -m $Best_fit -bb 1000 -redo -alrt 1000 -nt AUTO -pre iqtree


##########################################################
#plot
##########################################################

1. fastani_heatmap.R
2. iqtree_crc.R
3. iqtree_GMGC.R
4. gene_presence_absence_PCA.R




