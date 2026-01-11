## script for analyzing the Mutect gene variants


# filter = PASS
bcftools view -f PASS Z4_vs_D19.filtered.vcf.gz -o Z4_vs_D19.filtered.passed.vcf.gz

# filter AF >= 0.1
bcftools view -i 'FMT/AF>=0.1' Z5_vs_D19.filtered.passed.vcf.gz -Oz -o Z5_vs_D19.filtered.passed.AF10.vcf.gz 
## NOTE: This command will filter both WT D19 and sample, please double check ** We only want when sample AF >= 0.1

#index vcf file
bcftools index Z5_vs_D19.filtered.passed.AF10.vcf.gz




### DO THIS ONCE ###
# for gff3 gene annotation file, convert to bed file format
## 
awk '$3=="gene" {
  OFS="\t";
  # $1 = chrom, $4-1 = start (0-based), $5 = end
  # $9 = attributes (ID, Name, etc.)
  # column 5 = score (0), column 6 = strand ($7)
  print $1, $4-1, $5, $9, 0, $7
}' ../../Reference/FoxMRL8996_GeneModels_FilteredModels1.gff3 > FoxMRL8996_genes.strand.bed

# sort
bedtools sort \
  -i FoxMRL8996_genes.strand.bed \
  > FoxMRL8996_genes.strand.sorted.bed

### FINISH DO THIS ONCE ###


## For each individual file

# Convert VCF to BED (0-based)
bcftools query -f '%CHROM\t%POS0\t%POS\t%REF\t%ALT\t%FILTER\n' Z4_vs_D19.filtered.passed.AF10.vcf.gz > variant_to_gene/Z4_vs_D19.AF10.bed

cd variant_to_gene

# add a dummy strand info
awk 'BEGIN{OFS="\t"} {print $0, "."}' Z4_vs_D19.AF10.bed \
> Z4_vs_D19.AF10.strand.bed

# sort based on chromosome number
bedtools sort \
	-i Z4_vs_D19.AF10.strand.bed \
	> Z4_vs_D19.AF10.strand.sorted.bed


## find closest gene and strand info, up/downstream and other info
bedtools closest \
  -a Z4_vs_D19.AF10.strand.sorted.bed \
  -b FoxMRL8996_genes.strand.sorted.bed \
  -D a \
  > Z4_vs_D19.AF10.closest_gene.txt
  
  ## annotation this file ## within 1kb is considered ## doesn't distinguish 5'UTR, 3'UTR, Intron region
  awk 'BEGIN{
    OFS="\t"; 
    print "chr","start","end","ref","alt","FILTER","var_score","gene_chr","gene_start","gene_end","gene_info","gene_score","gene_strand","distance","location"
}
{
    dist = $14
    abs_dist = (dist < 0) ? -dist : dist  # absolute value

    if(dist == 0) {
        loc="CDS"
    } else if(abs_dist <= 1000) {
        loc = (dist < 0) ? "UPSTREAM" : "DOWNSTREAM"
    } else {
        loc="INTERGENIC"
    }
    
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,loc
}' Z4_vs_D19.AF10.closest_gene.txt > Z4_vs_D19.AF10.closest_gene.annotated_v2.txt



