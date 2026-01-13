# Load module on Unity
module load snpeff/2017-11-24

# modify config file
cp $SNPEFF_SPACK_ROOT/bin/snpEff.config snpeff.config.tmp
echo "FoxMRL8996.genome : FoxMRL8996" >> snpeff.config.tmp


# =============================
# Build the database
# =============================

echo "Building SnpEff database..."
## after moving fasta file and gff3 file into data/FoxMRL8996 and rename them
snpEff build -gff3 -v -c snpeff.config.tmp FoxMRL8996

# =============================
# Annotate VCF
# =============================

snpEff -v -c snpeff.config.tmp FoxMRL8996 Z5_vs_D19.filtered.passed.AF10.vcf.gz > Z5_vs_D19.filtered.passed.AF10.snpeff.vcf

mv snpEff_genes.txt outputs/Z5_snpEff_genes.txt
mv snpEff_summary.html outputs/Z5_snpEff_summary.html