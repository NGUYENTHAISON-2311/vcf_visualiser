# VCF Visualiser v2.1 - Complete README

[![Python](https://img.shields.io/badge/Python-3.6%2B-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](https://github.com)

> **Production-ready pipeline for variant annotation and visualization**  
> From VCF to publication-quality IGV heatmaps in one command

---

## 📋 Table of Contents

- [Overview](#overview)
- [Quick Start](#quick-start)
- [Installation](#installation)
- [Usage](#usage)
  - [CLI Mode](#cli-mode)
  - [Config Mode](#config-mode)
- [Input Requirements](#input-requirements)
- [Output Structure](#output-structure)
- [Test Datasets](#test-datasets)
- [Examples](#examples)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)
- [License](#license)

---

## 🎯 Overview

**VCF Visualiser** is an integrated pipeline that transforms VCF files into publication-ready visualizations for IGV (Integrative Genomics Viewer). It provides:

- **Comprehensive annotation** - Gene-based variant annotation
- **Aesthetic BED files** - Rich, informative variant annotations
- **Automated heatmaps** - Red-gradient density visualizations
- **BAM integration** - Optional coverage and allele frequency extraction
- **Full reproducibility** - Complete logging with MD5 checksums

### What Makes It Different?

✅ **No database required** - Works with any GFF file
✅ **100% offline** - No internet dependency
✅ **Single script** - No installation needed
✅ **Two simple modes** - CLI for quick runs, Config for projects
✅ **IGV-ready outputs** - Load directly in genome browser

---

## 🚀 Quick Start

### 30-Second Demo

```bash
# Download
wget https://raw.githubusercontent.com/NGUYENTHAISON-2311/vcf_visualiser/main/vcf_visualiser2.py
chmod +x vcf_visualiser2.py

# Run
python vcf_visualiser2.py \
  -i sample.vcf \
  -g genes.gff \
  -o results

# Load in IGV
# 1. File → Load from File → results/bed_files/sample.bed
# 2. File → Load from File → results/heatmaps/*.bedGraph
# 3. Right-click heatmap → Set Display Mode → Heatmap
```

---

## 📦 Installation

### Requirements

**Required:**
- Python 3.6 or higher (standard library only)

**Optional:**
- `samtools` - For BAM coverage/AF extraction

### Install

**Option 1: Direct Download**
```bash
wget https://raw.githubusercontent.com/NGUYENTHAISON-2311/vcf_visualiser/main/vcf_visualiser2.py
chmod +x vcf_visualiser2.py
```

**Option 2: Git Clone**
```bash
git clone https://github.com/NGUYENTHAISON-2311/vcf_visualiser.git
cd vcf_visualiser
```

### Optional: Install samtools

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install samtools
```

**macOS:**
```bash
brew install samtools
```

**Windows:**
```powershell
# Download from: https://github.com/samtools/samtools/releases
# Extract and add to PATH
```

**Conda (all platforms):**
```bash
conda install -c bioconda samtools
```

### Verify Installation

```bash
python vcf_visualiser2.py --version
# Output: v2.1.0

samtools --version  # Optional
# Output: samtools 1.x
```

---

## 💻 Usage

VCF Visualiser has **two modes**: CLI (simple) and Config (advanced).

### CLI Mode

**Use when:**
- Quick one-off analysis
- Single file or folder
- BAM files match VCF names

**Syntax:**
```bash
python vcf_visualiser2.py -i INPUT -g GFF -o OUTPUT [OPTIONS]
```

**Required Arguments:**
- `-i, --input` - VCF file or folder
- `-g, --gff` - GFF/GFF3/GTF annotation file
- `-o, --output` - Output directory

**Optional Arguments:**
- `-q, --min-quality` - Minimum QUAL score (default: 0)
- `-c, --min-coverage` - Minimum BAM coverage (default: 0)
- `-m, --color-mode` - BED color scheme: position, type, coverage, allele_freq (default: position)
- `-w, --window-size` - Heatmap window size in bp (default: 500)

#### CLI Examples

**Single file:**
```bash
python vcf_visualiser2.py \
  -i data/tumor.vcf \
  -g data/genes.gff \
  -o results
```

**Folder (auto-detect BAMs):**
```bash
python vcf_visualiser2.py \
  -i vcf_folder/ \
  -g genes.gff \
  -o results
```

**With filtering:**
```bash
python vcf_visualiser2.py \
  -i data/variants.vcf \
  -g data/genes.gff \
  -o results \
  -q 30 -c 50
```

**Custom color mode:**
```bash
python vcf_visualiser2.py \
  -i tumor.vcf \
  -g genes.gff \
  -o results \
  -m coverage \
  -w 1000
```

---

### Config Mode

**Use when:**
- Multi-sample projects
- BAM files have different names
- Need reproducible configuration
- Batch processing

**Syntax:**
```bash
python vcf_visualiser2.py --config CONFIG.json [OPTIONS]
```

#### Generate Config Template

```bash
python vcf_visualiser2.py --generate-config
```

This creates `vcf_visualiser_config.json`:

```json
{
  "gff_file": "path/to/annotations.gff",
  "output": "path/to/output",
  "samples": [
    {
      "name": "sample1",
      "vcf": "path/to/sample1.vcf",
      "bam": "path/to/sample1_reads.bam"
    },
    {
      "name": "sample2",
      "vcf": "path/to/sample2.vcf",
      "bam": "path/to/sample2_alignment.bam"
    },
    {
      "name": "sample3_no_bam",
      "vcf": "path/to/sample3.vcf"
    }
  ],
  "min_quality": 20,
  "min_coverage": 30,
  "color_mode": "position",
  "window_size": 500
}
```

#### Config Examples

**Basic config:**
```json
{
  "gff_file": "genes.gff",
  "output": "results",
  "samples": [
    {"vcf": "sample1.vcf"},
    {"vcf": "sample2.vcf"}
  ]
}
```

**With custom BAM names:**
```json
{
  "gff_file": "genes.gff",
  "output": "tumor_analysis",
  "samples": [
    {
      "name": "patient1_tumor",
      "vcf": "vcfs/p1_tumor.vcf",
      "bam": "bams/patient1_tumor_hg38_sorted.bam"
    },
    {
      "name": "patient1_normal",
      "vcf": "vcfs/p1_normal.vcf",
      "bam": "bams/patient1_normal_hg38_sorted.bam"
    }
  ],
  "min_quality": 30,
  "min_coverage": 50
}
```

**Run with config:**
```bash
python vcf_visualiser2.py --config project.json
```

**Override config with CLI:**
```bash
# Config has min_quality=20, override to 30
python vcf_visualiser2.py --config project.json -q 30

# Shows: CLI OVERRIDES APPLIED: min_quality: 30.0
```

---

## 📁 Input Requirements

### Required Files

#### 1. VCF File(s)

**Format:** Standard VCF format (v4.x)

```
##fileformat=VCFv4.2
#CHROM  POS     ID      REF  ALT  QUAL  FILTER  INFO
chr1    1000    .       A    G    45    PASS    DP=100;AF=0.45
chr1    2000    rs123   C    T    60    PASS    DP=80;AF=0.52
```

**Requirements:**
- Must have standard VCF columns
- Uncompressed (`.vcf` not `.vcf.gz`)
- Optional: INFO fields for DP, AF

#### 2. GFF Annotation File

**Formats:** GFF, GFF3, or GTF

```
##gff-version 3
chr1  source  gene  1000  5000  .  +  .  ID=gene1;Name=GENE1
chr1  source  CDS   1500  2500  .  +  .  Parent=gene1
```

**Requirements:**
- Must have gene features
- Name/ID attributes for gene names

#### 3. BAM Files (Optional)

**Format:** Aligned BAM files

**Requirements:**
- Same reference as VCF
- **CLI mode:** BAM name must match VCF name
  - `sample.vcf` → `sample.bam` ✅
  - `sample.vcf` → `sample_reads.bam` ❌
- **Config mode:** Any BAM name (explicit pairing)
- Index (`.bai`) not required

### File Structure Examples

**CLI Mode:**
```
project/
├── vcf_data/
│   ├── sample1.vcf
│   ├── sample1.bam         ← Matches!
│   ├── sample2.vcf
│   ├── sample2.bam         ← Matches!
│   └── genes.gff
```

**Config Mode:**
```
project/
├── vcfs/
│   ├── patient1_tumor.vcf
│   └── patient2_tumor.vcf
├── bams/
│   ├── p1_tumor_aligned.bam    ← Different name OK!
│   └── p2_tumor_aligned.bam    ← Different name OK!
└── reference/
    └── genes.gff
```

---

## 📤 Output Structure

```
results/
├── bed_files/                    # Load in IGV
│   ├── sample1.bed
│   └── sample2.bed
│
├── tables/                       # TSV for analysis
│   ├── sample1_details.tsv
│   └── sample2_details.tsv
│
├── summaries/                    # Statistics
│   ├── sample1_summary.txt
│   ├── sample2_summary.txt
│   └── combined_summary.txt
│
├── heatmaps/                     # Load in IGV
│   ├── sample1_gene_hotspots.bedGraph
│   ├── sample1_density_500bp.bedGraph
│   ├── sample1_af_landscape.bedGraph   # If BAM available
│   ├── sample2_gene_hotspots.bedGraph
│   └── ...
│
└── processing_log.json           # Reproducibility
```

### Output Files Explained

#### BED Files (`bed_files/`)

**Purpose:** Color-coded variants for IGV

**Format:** BED9 with RGB colors

**Annotations include:**
- Variant (ref>alt)
- Type (SNV/INS/DEL)
- Gene name(s)
- Location (coding/exon/gene/intergenic)
- Quality score
- Coverage (if BAM)
- Allele frequencies
- Genomic position

**Example:**
```
chr1  999  1000  Var:A>G | Type:SNV | Gene:BRCA1 | Loc:coding | QUAL:85 | Cov:150x | BAM_AF:0.93 | Pos:chr1:1000  1000  +  999  1000  255,0,0
```

#### TSV Tables (`tables/`)

**Purpose:** Structured data for analysis

**Columns:**
- Chrom, Pos, Ref, Alt, Type
- QUAL, Filter
- Gene, Location
- VCF_DP, VCF_AF
- BAM_Cov, BAM_AF, BAM_Ref, BAM_Alt

**Example:**
```
Chrom  Pos   Ref  Alt  Type  QUAL  Filter  Gene   Location  VCF_DP  VCF_AF  BAM_Cov  BAM_AF
chr1   1000  A    G    SNV   85    PASS    BRCA1  coding    145     0.950   150      0.933
```

#### Heatmaps (`heatmaps/`)

**Three types:**

1. **Gene Hotspots** (`*_gene_hotspots.bedGraph`)
   - Variant count per gene
   - Identifies frequently mutated genes
   - Red intensity = variant count

2. **Window Density** (`*_density_500bp.bedGraph`)
   - Variants per genomic window
   - Identifies clustering regions
   - Red intensity = variant density

3. **AF Landscape** (`*_af_landscape.bedGraph`)
   - Average allele frequency per window
   - Identifies LOH regions
   - Red intensity = AF (0-100%)
   - Only generated with BAM data

**All use RED color scheme for publication visibility**

#### Summaries (`summaries/`)

**Per-sample summary:**
```
======================================================================
SAMPLE: patient1
======================================================================

FILTERING:
  Min QUAL: 20
  Min Coverage: 30

VARIANTS:
  Total: 1247
  Passed: 856 (68.6%)

BY TYPE:
  SNV: 723
  INS: 89
  DEL: 44

BY LOCATION:
  coding: 234
  exon: 145
  gene: 312
  intergenic: 165
```

**Combined summary:**
```
======================================================================
COMBINED ANALYSIS SUMMARY
======================================================================

Samples: 3

PER-SAMPLE STATS:
----------------------------------------------------------------------
Sample              Total    Passed    Coding
----------------------------------------------------------------------
patient1             1247       856       234
patient2             1543       982       298
patient3              987       645       187
```

#### Processing Log (`processing_log.json`)

**Purpose:** Complete reproducibility record

```json
{
  "version": "2.1.0",
  "timestamp": "2024-02-15T10:30:00",
  "command": "python vcf_visualiser2.py -i vcf -g genes.gff -o results",
  "parameters": {
    "min_quality": "20",
    "min_coverage": "30",
    "color_mode": "position"
  },
  "input_files": {
    "vcf_sample1": {
      "path": "/full/path/sample1.vcf",
      "size": 12345,
      "md5": "abc123..."
    },
    "gff": {
      "path": "/full/path/genes.gff",
      "size": 54321,
      "md5": "def456..."
    }
  },
  "steps": [
    {"step": "File detection", "time": "..."},
    {"step": "Variant processing", "time": "..."}
  ]
}
```

---

## 🧪 Test Datasets

### Test Dataset 1: Basic SNP Analysis

**Description:** Small bacterial genome with 50 SNPs

**Download:**
```bash
wget https://github.com/your-repo/test-data/test1_basic.tar.gz
tar -xzf test1_basic.tar.gz
```

**Contents:**
```
test1_basic/
├── ecoli_variants.vcf      # 50 SNPs
├── ecoli_genome.gff        # E. coli genome annotation
└── README.txt
```

**Run:**
```bash
python vcf_visualiser2.py \
  -i test1_basic/ecoli_variants.vcf \
  -g test1_basic/ecoli_genome.gff \
  -o test1_output
```

**Expected output:**
- 50 total variants
- ~40 in genes
- Gene hotspots showing ~5 genes with multiple hits
- Runtime: ~5 seconds

---

### Test Dataset 2: Multi-Sample Cohort

**Description:** 5 samples from small genome (100 variants each)

**Download:**
```bash
wget https://github.com/your-repo/test-data/test2_cohort.tar.gz
tar -xzf test2_cohort.tar.gz
```

**Contents:**
```
test2_cohort/
├── vcfs/
│   ├── sample1.vcf
│   ├── sample2.vcf
│   ├── sample3.vcf
│   ├── sample4.vcf
│   └── sample5.vcf
├── genome.gff
└── config.json              # Pre-configured
```

**Run (CLI):**
```bash
python vcf_visualiser2.py \
  -i test2_cohort/vcfs/ \
  -g test2_cohort/genome.gff \
  -o test2_output
```

**Run (Config):**
```bash
python vcf_visualiser2.py --config test2_cohort/config.json
```

**Expected output:**
- 5 samples processed
- ~500 total variants
- Comparative heatmaps showing sample differences
- Runtime: ~30 seconds

---

### Test Dataset 3: With BAM Coverage

**Description:** 1 sample with VCF + BAM (realistic coverage data)

**Download:**
```bash
wget https://github.com/your-repo/test-data/test3_bam.tar.gz
tar -xzf test3_bam.tar.gz
```

**Contents:**
```
test3_bam/
├── tumor.vcf                # 200 variants
├── tumor.bam                # Aligned reads
├── genes.gff
└── README.txt
```

**Prerequisites:**
- samtools must be installed

**Run:**
```bash
python vcf_visualiser2.py \
  -i test3_bam/tumor.vcf \
  -g test3_bam/genes.gff \
  -o test3_output \
  -c 30                      # Use coverage filter
```

**Expected output:**
- BAM coverage extracted
- AF calculated from BAM
- AF landscape heatmap generated
- Coverage-based filtering applied
- Runtime: ~2 minutes

---

### Test Dataset 4: Human Chromosome Subset

**Description:** Human chr22 region with real annotations

**Download:**
```bash
wget https://github.com/your-repo/test-data/test4_human.tar.gz
tar -xzf test4_human.tar.gz
```

**Contents:**
```
test4_human/
├── chr22_variants.vcf       # 1000 variants
├── chr22_genes.gff          # GENCODE annotations
└── README.txt
```

**Run:**
```bash
python vcf_visualiser2.py \
  -i test4_human/chr22_variants.vcf \
  -g test4_human/chr22_genes.gff \
  -o test4_output \
  -q 30 -w 1000
```

**Expected output:**
- ~800 variants passing filter
- Gene hotspots in known cancer genes
- High-quality heatmaps
- Runtime: ~1 minute

---

### Test Dataset 5: Non-Model Organism

**Description:** Custom organism with minimal annotation

**Download:**
```bash
wget https://github.com/your-repo/test-data/test5_custom.tar.gz
tar -xzf test5_custom.tar.gz
```

**Contents:**
```
test5_custom/
├── variants.vcf             # 150 variants
├── custom_genome.gff        # Minimal annotation
├── sample.bam               # Optional
└── README.txt
```

**Run:**
```bash
python vcf_visualiser2.py \
  -i test5_custom/variants.vcf \
  -g test5_custom/custom_genome.gff \
  -o test5_output
```

**Purpose:** Demonstrates working with custom organisms

**Expected output:**
- Works without pre-built databases
- Annotations based on minimal GFF
- Runtime: ~10 seconds

---

### Test Dataset 6: Config Mode Advanced

**Description:** Complex project with different BAM names

**Download:**
```bash
wget https://github.com/your-repo/test-data/test6_config.tar.gz
tar -xzf test6_config.tar.gz
```

**Contents:**
```
test6_config/
├── vcfs/
│   ├── patient1_tumor.vcf
│   └── patient1_normal.vcf
├── alignments/
│   ├── p1_T_hg38_sorted.bam     # Different names!
│   └── p1_N_hg38_sorted.bam
├── reference/
│   └── genes.gff
└── config.json
```

**Config (`config.json`):**
```json
{
  "gff_file": "reference/genes.gff",
  "output": "paired_analysis",
  "samples": [
    {
      "name": "patient1_tumor",
      "vcf": "vcfs/patient1_tumor.vcf",
      "bam": "alignments/p1_T_hg38_sorted.bam"
    },
    {
      "name": "patient1_normal",
      "vcf": "vcfs/patient1_normal.vcf",
      "bam": "alignments/p1_N_hg38_sorted.bam"
    }
  ],
  "min_quality": 30,
  "min_coverage": 50,
  "color_mode": "allele_freq"
}
```

**Run:**
```bash
python vcf_visualiser2.py --config test6_config/config.json
```

**Purpose:** Demonstrates tumor-normal pairing with config mode

**Expected output:**
- Paired samples processed
- BAM coverage and AF extracted
- Side-by-side comparison possible in IGV
- Runtime: ~3 minutes

---

### Generate Your Own Test Data

**Quick test generator:**
```bash
cat > generate_test.py << 'EOF'
#!/usr/bin/env python3
import random

# Generate minimal VCF
with open('test.vcf', 'w') as f:
    f.write('##fileformat=VCFv4.2\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for i in range(100):
        pos = random.randint(1000, 10000)
        ref = random.choice('ACGT')
        alt = random.choice('ACGT')
        qual = random.randint(20, 100)
        f.write(f'chr1\t{pos}\t.\t{ref}\t{alt}\t{qual}\tPASS\tDP=50\n')

# Generate minimal GFF
with open('test.gff', 'w') as f:
    f.write('##gff-version 3\n')
    for i in range(10):
        start = i * 1000
        end = start + 500
        f.write(f'chr1\t.\tgene\t{start}\t{end}\t.\t+\t.\tID=gene{i};Name=GENE{i}\n')

print("Generated test.vcf and test.gff")
EOF

python3 generate_test.py
```

**Test:**
```bash
python vcf_visualiser2.py -i test.vcf -g test.gff -o test_output
```

---

## 📚 Examples

### Example 1: Single Sample Quick Analysis

```bash
# You have one VCF, want quick visualization
python vcf_visualiser2.py \
  -i sample.vcf \
  -g genes.gff \
  -o results

# Load in IGV:
# - results/bed_files/sample.bed
# - results/heatmaps/sample_gene_hotspots.bedGraph
# - results/heatmaps/sample_density_500bp.bedGraph
```

---

### Example 2: Cohort with Matching BAMs

```bash
# Folder structure:
# cohort/
#   ├── patient1.vcf
#   ├── patient1.bam
#   ├── patient2.vcf
#   ├── patient2.bam
#   └── genes.gff

python vcf_visualiser2.py \
  -i cohort/ \
  -g cohort/genes.gff \
  -o cohort_results \
  -q 20 -c 30

# Processes all samples with coverage filtering
```

---

### Example 3: Tumor-Normal Pair (Config Mode)

```bash
# Create config
cat > tumor_normal.json << 'EOF'
{
  "gff_file": "genes.gff",
  "output": "paired_analysis",
  "samples": [
    {
      "name": "tumor",
      "vcf": "tumor.vcf",
      "bam": "tumor_aligned.bam"
    },
    {
      "name": "normal",
      "vcf": "normal.vcf",
      "bam": "normal_aligned.bam"
    }
  ],
  "min_quality": 30,
  "min_coverage": 50,
  "color_mode": "allele_freq",
  "window_size": 1000
}
EOF

# Run
python vcf_visualiser2.py --config tumor_normal.json

# Compare in IGV side-by-side
```

---

### Example 4: Filtering Workflow

```bash
# Step 1: Pre-filter with bcftools
bcftools view -i 'QUAL>20 && DP>10' raw.vcf > filtered.vcf

# Step 2: Visualize
python vcf_visualiser2.py \
  -i filtered.vcf \
  -g genes.gff \
  -o filtered_viz

# Step 3: Review in IGV
# Load filtered_viz/bed_files/*.bed
```

---

### Example 5: Non-Model Organism

```bash
# Custom genome - no databases available
python vcf_visualiser2.py \
  -i custom_species.vcf \
  -g custom_annotation.gff \
  -o custom_results

# Works immediately - no database download!
```

---

### Example 6: Large Project with Overrides

```bash
# Use config for base parameters
cat > project.json << 'EOF'
{
  "gff_file": "genes.gff",
  "output": "default_output",
  "samples": [ ... ],
  "min_quality": 20
}
EOF

# Try different thresholds
python vcf_visualiser2.py --config project.json -q 10 -o permissive/
python vcf_visualiser2.py --config project.json -q 30 -o stringent/
python vcf_visualiser2.py --config project.json -q 50 -o very_stringent/

# Compare results
```

---

## 🔧 Troubleshooting

### Issue: "samtools not found"

**Symptom:**
```
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WARNING: samtools not found in PATH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
```

**Impact:**
- BAM coverage NOT calculated
- BAM_AF column shows N/A
- AF landscape may fail (unless VCF has AF)

**Solution:**
```bash
# Install samtools
sudo apt-get install samtools  # Ubuntu
brew install samtools           # macOS
conda install -c bioconda samtools  # Conda

# Verify
samtools --version
```

**Workaround:**
- Use VCF_AF if available
- Skip BAM features

---

### Issue: "Coverage filter requires BAM files"

**Symptom:**
```
error: Coverage filter (-c) requires BAM files
```

**Cause:**
Config/input has no BAM files but `-c` flag used

**Solution:**
```bash
# Remove -c flag if no BAM
python vcf_visualiser2.py -i data.vcf -g genes.gff -o results

# OR add BAM files
python vcf_visualiser2.py -i data.vcf -g genes.gff -o results -c 30
# (Ensure data.bam exists)
```

---

### Issue: "No VCF files found"

**Symptom:**
```
error: No VCF files found in VCF directory
```

**Cause:**
- Wrong directory
- Files are `.vcf.gz` (compressed)

**Solution:**
```bash
# Check directory
ls vcf_folder/*.vcf

# Decompress if needed
gunzip vcf_folder/*.vcf.gz

# Verify
ls vcf_folder/*.vcf
```

---

### Issue: "GFF file not found"

**Symptom:**
```
error: GFF file not found: genes.gff
```

**Solution:**
```bash
# Use absolute path
python vcf_visualiser2.py \
  -i data.vcf \
  -g /full/path/to/genes.gff \
  -o results

# Or check file exists
ls -l genes.gff
```

---

### Issue: BAM detected but no AF calculated

**Symptom:**
- BAM shown in scan
- AF landscape fails: "✗ (no AF data)"
- TSV BAM_AF column is N/A

**Cause:**
samtools not available

**Solution:**
Install samtools (see above)

**Check:**
```bash
# View TSV
head results/tables/sample_details.tsv
# Check BAM_AF column - should have numbers not N/A
```

---

### Issue: "Heatmap generation failed"

**Symptom:**
```
[sample]
  Gene hotspots... ✗
```

**Possible causes:**
1. No variants in GFF regions
2. GFF parsing error
3. Chromosome name mismatch

**Debug:**
```bash
# Check VCF chromosomes
grep -v "^#" sample.vcf | cut -f1 | sort -u

# Check GFF chromosomes  
grep -v "^#" genes.gff | cut -f1 | sort -u

# Must match!
```

---

### Issue: Out of memory

**Symptom:**
```
Killed
```

**Cause:**
Too many variants or large BAM files

**Solution:**
```bash
# Filter first
bcftools view -i 'QUAL>30' large.vcf > filtered.vcf
python vcf_visualiser2.py -i filtered.vcf -g genes.gff -o results

# Or process in batches
bcftools view -r chr1 large.vcf > chr1.vcf
python vcf_visualiser2.py -i chr1.vcf -g genes.gff -o chr1_results
```

---

## 📖 Citation

If you use VCF Visualiser in your research, please cite:

```bibtex
@software{vcf_visualiser,
  title = {VCF Visualiser: Integrated Pipeline for Variant Annotation and Visualization},
  author = {Thai Son NGUYEN},
  year = {2026},
  version = {2.1.0},
  url = {https://github.com/your-repo/vcf-visualizer}
}
```

---

## 📄 License

MIT License - see LICENSE file for details

---

## 🤝 Contributing

Contributions welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

---


**Made with 🧬 for the genomics community**
