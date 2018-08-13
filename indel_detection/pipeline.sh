#restrict to genes which are at least in 50 isolate genomes
awk -v threshold="50" -f ~/pseudo_genomics/src/indel_detection/filter.awk  < gene_presence_absence.csv  > core_genes_50.txt
core_genes=~/pseudo_genomics/results/assembly/v2/roary/v5/out_95/core_genes_50.txt
#number of core genes
core_genes_wc=$(wc -l $core_genes | cut -f1 -d" ")
#extract fasta sequences for each gene family
python  ~/pseudo_genomics/src/PseudoGenomics/indel_detection/gene_clusters2multi_fasta.py extracted_proteins_nt/ <(head -n $core_genes_wc clustered_proteins)  ~/pseudo_genomics/results/assembly/v2/prokka/v4/*/*.ffn
#cd to roary dir 
mkdir extracted_protein_nt
cd extracted_protein_nt
#multiple sequence alignment
for i in `ls`; do i=`echo $i | cut -f1 -d "."` ; echo "mafft $i.fasta > $i.aln.fasta"; done | parallel -j 20
cd ..
mkdir indels
cd indels
#variant calling on msa
cat $core_genes | while read i; do echo "java -jar ~/software/jvarkit/dist/msa2vcf.jar < ~/pseudo_genomics/results/assembly/v2/roary/v5/out_95/extracted_proteins_nt/${i}.aln.fasta > ${i}.vcf"; done
#vcf to indel yes/no vector, stats and gff
cat $core_genes | while read i ; do echo "python /net/sgi/metagenomics/projects/pseudo_genomics/src/PseudoGenomics/indel_detection/vcf2indel.py $i.vcf $i ${i}_indels.txt ${i}_indels.gff ${i}_indel_stats.txt"; done

python /net/sgi/metagenomics/projects/pseudo_genomics/src/PseudoGenomics/indel_detection/generate_indel_features.py 
<(cut -f1 ~/pseudo_genomics/results/assembly/v2/roary/v5/out_95/gene_presence_absence.Rtab | tail -n+2 | grep -v hdl )   
~/pseudo_genomics/results/classification/v6/gpa/out_95/annot.txt 
indel_annot.txt 
indels_stats.txt 
~/pseudo_genomics/results/assembly/v2/abricate2roary/ncbi/roary_PA14_abricate.txt


1. target genes
2. gene families
3. alignment
4. vcf
5. only indels
6. indels table
