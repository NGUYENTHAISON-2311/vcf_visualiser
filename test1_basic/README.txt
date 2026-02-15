TEST DATASET 1: Basic SNP Analysis
===================================

Description:
- Small E. coli genome subset
- 50 SNPs across 9 genes
- Simple test case for learning

Files:
- ecoli_variants.vcf   (50 SNPs)
- ecoli_genome.gff     (9 genes)

Run:
python vcf_visualizer_v2.1.py \
  -i ecoli_variants.vcf \
  -g ecoli_genome.gff \
  -o test1_output

Expected:
- 50 variants total
- All in genes
- ~5-6 variants per gene
- Gene hotspots showing lacZ, rpoB, recA with multiple hits
- Runtime: ~5 seconds

Good for:
- First-time users
- Testing installation
- Learning the tool
