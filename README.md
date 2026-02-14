# VCF Visualizer v1.0 - Complete Documentation

## 🎯 Overview

**VCF Visualizer** is a comprehensive, production-ready tool for generating annotated BED files from VCF, BAM, and GFF files for visualization in IGV.

### Key Features

✅ **Multiple VCF Processing** - Process one or many VCF files  
✅ **Contig Validation** - Checks chromosome name compatibility  
✅ **Multiple Color Modes** - Position, type, allele frequency, quality  
✅ **Flexible Filtering** - Quality and coverage filters  
✅ **Comprehensive Output** - BED, summaries, tables  
✅ **Full Reproducibility** - Complete logging with checksums  
✅ **Production Ready** - Error handling, validation, documentation  

## 📁 Input/Output Structure

### Input

```
input_folder/
├── sample1.vcf      ← One or more VCF files
├── sample2.vcf
├── sample3.vcf
├── annotations.gff  ← Single GFF file (required)
└── sample1.bam      ← Optional BAM files (matched by name)
```

### Output

```
output_folder/
├── bed_files/
│   ├── sample1.bed          ← Load these in IGV!
│   ├── sample2.bed
│   └── sample3.bed
├── summaries/
│   ├── sample1_summary.txt  ← Human-readable summaries
│   ├── sample2_summary.txt
│   └── sample3_summary.txt
├── tables/
│   ├── sample1_details.tsv  ← Detailed variant tables
│   ├── sample2_details.tsv
│   └── sample3_details.tsv
└── processing_log.json      ← Complete processing log
```

## 🎨 Color Modes

### 1. Position Mode (Default)

Colors based on genomic location - **works without BAM**

| Color | Location | Priority |
|-------|----------|----------|
| 🔴 Red | Coding (CDS) | Highest |
| 🟠 Orange | Exon | High |
| 🔵 Blue | Gene (intronic) | Medium |
| 🟢 Green | Other genic | Low |
| ⚫ Gray | Intergenic | Lowest |

```bash
python vcf_visualizer.py -i data -o results -m position
```

### 2. Type Mode

Colors based on variant type - **works without BAM**

| Color | Type | Description |
|-------|------|-------------|
| 🔵 Blue | SNV | Single nucleotide variant |
| 🟢 Green | INS | Insertion |
| 🔴 Red | DEL | Deletion |
| 🟣 Magenta | MNV | Multi-nucleotide variant |

```bash
python vcf_visualizer.py -i data -o results -m type
```

### 3. Allele Frequency Mode

Colors based on allele frequency - **requires -b BAM**

| Color | AF Range | Interpretation |
|-------|----------|----------------|
| ⚪ Light gray | 0-0.1 | Very rare |
| 🔵 Light blue | 0.1-0.25 | Rare |
| 🔵 Blue | 0.25-0.4 | Low frequency |
| 🟢 Green | 0.4-0.6 | Medium (heterozygous) |
| 🔴 Red | 0.6-1.0 | High (homozygous) |

```bash
python vcf_visualizer.py -i data -o results -b sample.bam -m allele_freq
```

### 4. Quality Mode

Colors based on quality evidence - **requires -b BAM**

| Color | QUAL Range | Confidence |
|-------|-----------|------------|
| 🔴 Light red | <20 | Low |
| 🟡 Yellow | 20-40 | Medium |
| 🟢 Light green | 40-60 | High |
| 🟢 Dark green | >60 | Very high |

```bash
python vcf_visualizer.py -i data -o results -b sample.bam -m quality
```

## 💻 Command-Line Options

### Required

```
-i, --input-dir DIR       Input folder with VCF, GFF, (BAM) files
-o, --output-dir DIR      Output folder
```

### Optional - Filtering

```
-q, --min-quality FLOAT   Minimum QUAL score (default: 0)
-c, --min-coverage INT    Minimum coverage (requires -b, default: 0)
```

### Optional - BAM Processing

```
-b, --bam FILE            Process single VCF with this BAM
                          Enables allele_freq and quality modes
```

### Optional - Color Mode

```
-m, --color-mode MODE     Color mode (default: position)
                          Choices: position, type, allele_freq, quality
```

### Optional - Reproducibility

```
-r, --reference NAME      Reference genome name (for documentation)
--check-contigs           Verify contig name compatibility
--version                 Show version number
--save-params FILE        Save parameters to JSON file
--load-params FILE        Load parameters from JSON file
```

## 📋 Usage Examples

### Example 1: Default - Multiple VCFs, Position Coloring

```bash
python vcf_visualizer.py \
  -i cohort_data \
  -o cohort_results \
  -q 20

# Processes all VCFs
# Colors by genomic position
# Filters QUAL >= 20
```

### Example 2: Single VCF with BAM, AF Coloring

```bash
python vcf_visualizer.py \
  -i patient_data \
  -o patient_results \
  -b patient_001.bam \
  -m allele_freq \
  -q 30 \
  -c 50

# Processes one VCF with matching BAM
# Colors by allele frequency
# Filters QUAL >= 30 AND Coverage >= 50
```

### Example 3: Variant Type Analysis

```bash
python vcf_visualizer.py \
  -i variants \
  -o results \
  -m type \
  --check-contigs \
  -r hg38

# Colors by variant type (SNV/INS/DEL)
# Checks contig compatibility
# Documents reference genome
```

### Example 4: High-Stringency Clinical

```bash
python vcf_visualizer.py \
  -i clinical_samples \
  -o clinical_results \
  -b patient.bam \
  -m quality \
  -q 50 \
  -c 100 \
  -r GRCh38

# Very strict quality filters
# Colors by quality evidence
# Clinical-grade filtering
```

## 📊 BED File Format

Each variant is annotated with comprehensive information:

```
chr17  7576999  7577000  Var:G>A | Type:SNV | Gene:TP53 | Loc:coding | QUAL:70.0 | BAM_AF:0.93 | VCF_AF:0.95 | Pos:chr17:7577000  350  -  7576999  7577000  255,0,0
```

**BED Fields:**
1. Chromosome
2. Start (0-based)
3. End
4. Name (comprehensive annotation)
5. Score (for sorting)
6. Strand
7. ThickStart
8. ThickEnd
9. RGB color

## 📈 Output Files Explained

### 1. BED Files (`bed_files/`)

**Purpose:** Load directly into IGV for visualization

**Format:** Standard BED12 with RGB colors

**Contains:**
- Variant change (Ref>Alt)
- Variant type (SNV/INS/DEL)
- Gene location
- Quality scores
- Allele frequencies
- Position

**Usage in IGV:**
```
File → Load from File → sample.bed
```

### 2. Summary Files (`summaries/`)

**Purpose:** Human-readable analysis summary

**Contains:**
- Sample information
- Filtering criteria
- Variant counts (total, passed, failed)
- Statistics by type, location, contig

**Example:**
```
Total variants: 247
Passed filters: 189 (76.5%)
Failed quality: 45
Failed coverage: 13

BY VARIANT TYPE:
SNV: 230 (93.1%)
INS: 12 (4.9%)
DEL: 5 (2.0%)

BY LOCATION:
coding: 95 (38.5%)
exon: 57 (23.1%)
gene: 25 (10.1%)
```

### 3. Detail Tables (`tables/`)

**Purpose:** Detailed variant information for downstream analysis

**Format:** Tab-separated values (TSV)

**Columns:**
- Chrom, Pos, Ref, Alt
- Type, QUAL, Filter
- Gene, Location
- VCF_DP, VCF_AF
- BAM_Cov, BAM_AF, BAM_Ref, BAM_Alt

**Use with:**
- Excel/LibreOffice
- R/Python data analysis
- Database import
- Custom scripts

### 4. Processing Log (`processing_log.json`)

**Purpose:** Complete reproducibility record

**Contains:**
- Version number
- Timestamp
- Complete command line
- All parameters
- Input file information:
  - Full paths
  - File sizes
  - MD5 checksums
  - Modification times
- Processing steps
- Statistics

**Benefits:**
- Verify input files haven't changed
- Reproduce exact analysis
- Track processing workflow
- Document parameters

## 🔍 Contig Compatibility Check

Use `--check-contigs` to verify chromosome names match:

```bash
python vcf_visualizer.py -i data -o results --check-contigs
```

**Output:**
```
VCF contigs: 25
GFF contigs: 25
Common contigs: 25
Overlap: 100.0%
✓ Compatible
```

**Common Issues:**
- VCF uses "chr1", GFF uses "1"
- Different reference builds
- Extra contigs in one file

**Solution:** Ensure VCF and GFF use same reference and naming convention

## 🔬 Processing Modes

### Multiple VCF Mode (Default)

**Trigger:** No -b option specified

**Behavior:**
- Processes all VCFs in input folder
- Uses position coloring only
- No BAM data
- Generates outputs for each VCF

**Best for:**
- Cohort studies
- Multiple samples
- Quick visualization
- Position-based analysis

### Single VCF + BAM Mode

**Trigger:** -b option specified

**Behavior:**
- Processes one VCF with matched BAM
- All color modes available
- Full BAM annotation
- Advanced filtering

**Best for:**
- Deep analysis of single sample
- Quality assessment
- Allele frequency analysis
- Clinical validation

## 📝 Real-World Workflows

### Workflow 1: Cohort Variant Discovery

```bash
# Input: 50 patient VCFs
cohort/
├── patient_001.vcf through patient_050.vcf
└── genes.gff

# Run
python vcf_visualizer.py \
  -i cohort \
  -o cohort_analysis \
  -q 20 \
  -m position \
  --check-contigs \
  -r hg38

# Output: 50 BED files, one per patient
# Load all in IGV to compare
```

### Workflow 2: Clinical Sample QC

```bash
# Input: Single patient with BAM
clinical/
├── patient_tumor.vcf
├── patient_tumor.bam
└── clinical_genes.gff

# Run
python vcf_visualizer.py \
  -i clinical \
  -o clinical_qc \
  -b patient_tumor.bam \
  -m quality \
  -q 50 \
  -c 100 \
  -r GRCh38

# Review quality-colored variants in IGV
```

### Workflow 3: Somatic Variant Analysis

```bash
# Input: Tumor sample
tumor/
├── somatic_variants.vcf
├── tumor_reads.bam
└── cancer_genes.gff

# Run
python vcf_visualizer.py \
  -i tumor \
  -o somatic_analysis \
  -b tumor_reads.bam \
  -m allele_freq \
  -q 20 \
  -c 50

# Identify low-frequency somatic mutations
```

## 🛠️ Reproducibility Features

### 1. Version Tracking

```bash
python vcf_visualizer.py --version
# VCF Visualizer v1.0.0
```

### 2. Complete Command Logging

All parameters saved in `processing_log.json`:
```json
{
  "command_line": "vcf_visualizer.py -i data -o results -q 20",
  "parameters": {
    "input_dir": "data",
    "output_dir": "results",
    "min_quality": "20.0"
  }
}
```

### 3. Input File Verification

MD5 checksums for all input files:
```json
{
  "input_files": {
    "vcf_sample": {
      "path": "/full/path/sample.vcf",
      "checksum_md5": "abc123...",
      "size_bytes": 12345
    }
  }
}
```

**Verify files haven't changed:**
```bash
md5sum sample.vcf
# Compare with checksum in log
```

### 4. Reference Documentation

```bash
python vcf_visualizer.py -i data -o results -r "GRCh38_v42"
```

Logged in `processing_log.json`:
```json
{
  "parameters": {
    "reference": "GRCh38_v42"
  }
}
```

## ⚙️ Advanced Options

### Save/Load Parameters

**Save parameters for reuse:**
```bash
python vcf_visualizer.py \
  -i data -o results \
  -q 30 -c 50 -m allele_freq \
  --save-params analysis_params.json
```

**Load saved parameters:**
```bash
python vcf_visualizer.py --load-params analysis_params.json
```

### Custom Filtering Strategies

**Permissive (exploratory):**
```bash
-q 10 -c 10
```

**Standard (balanced):**
```bash
-q 20 -c 30
```

**Stringent (clinical):**
```bash
-q 50 -c 100
```

## 🔧 Troubleshooting

### No BAM Coverage

**Issue:** BAM_AF shows N/A

**Cause:** samtools not installed or BAM not provided

**Solution:**
```bash
# Install samtools
sudo apt-get install samtools

# Or run without BAM
python vcf_visualizer.py -i data -o results
```

### Contig Mismatch

**Issue:** "No common contigs"

**Cause:** VCF and GFF use different chromosome names

**Solution:** Ensure same reference genome and naming

### Multiple GFF Files

**Issue:** "Multiple GFF files found"

**Solution:** Tool uses first one - remove extras or specify which to keep

## 📊 Performance

**Typical runtimes:**

| Variants | Without BAM | With BAM |
|----------|-------------|----------|
| 1,000 | ~1 sec | ~5 sec |
| 10,000 | ~5 sec | ~45 sec |
| 100,000 | ~45 sec | ~7 min |

**BAM processing is the bottleneck** - coverage extraction takes time

**Tip:** For large datasets without needing BAM AF, skip -b option

## 🎓 Best Practices

1. **Always check contigs first:**
   ```bash
   --check-contigs
   ```

2. **Document reference:**
   ```bash
   -r GRCh38
   ```

3. **Start permissive, then filter:**
   ```bash
   # First run
   -q 0 -c 0
   
   # Check summaries, then apply filters
   -q 20 -c 30
   ```

4. **Use appropriate color mode:**
   - Multiple samples → position
   - Single sample deep dive → allele_freq or quality

5. **Save processing logs** for reproducibility

## 📦 Dependencies

- Python 3.6+
- samtools (optional, for BAM processing)

**No Python package dependencies!** Uses only standard library.

## 🎯 Summary

VCF Visualizer provides:

✅ **Flexible processing** - One or many VCFs  
✅ **Multiple color modes** - Position, type, AF, quality  
✅ **Comprehensive output** - BED, summaries, tables  
✅ **Full reproducibility** - Logs, checksums, versions  
✅ **Production ready** - Error handling, validation  
✅ **Easy to use** - Clear options, good defaults  

**One tool. Complete solution. IGV ready.**