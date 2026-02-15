#!/usr/bin/env python3
"""
VCF Visualizer v2.1 - Complete Variant Visualization Pipeline
Version: 2.1.0 (Updated)
"""

import argparse
import sys
import json
import hashlib
import shutil
from pathlib import Path
from typing import List, Dict, Optional
from datetime import datetime
from collections import defaultdict
import subprocess
import re

__version__ = "2.1.0"

# ============================================================================
# UTILITY CLASSES
# ============================================================================

class Logger:
    """Reproducibility logger"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.log = {
            'version': __version__,
            'timestamp': datetime.now().isoformat(),
            'command': ' '.join(sys.argv),
            'parameters': {},
            'input_files': {},
            'steps': []
        }
    
    def log_param(self, key: str, value):
        self.log['parameters'][key] = str(value)
    
    def log_file(self, label: str, filepath: Path):
        if filepath and filepath.exists():
            self.log['input_files'][label] = {
                'path': str(filepath.absolute()),
                'size': filepath.stat().st_size,
                'md5': self._md5(filepath)
            }
    
    def log_step(self, step: str):
        self.log['steps'].append({'step': step, 'time': datetime.now().isoformat()})
    
    def _md5(self, filepath: Path) -> str:
        hash_md5 = hashlib.md5()
        with open(filepath, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    
    def save(self):
        with open(self.output_dir / "processing_log.json", 'w') as f:
            json.dump(self.log, f, indent=2)

class BAMParser:
    """Integrated BAM parser"""
    
    def __init__(self):
        self.samtools_available = self._check_samtools()
        if not self.samtools_available:
            self.warning_shown = False
    
    def _check_samtools(self) -> bool:
        try:
            subprocess.run(['samtools', '--version'], capture_output=True, timeout=2)
            return True
        except:
            return False
    
    def show_warning(self):
        """Show samtools warning once"""
        if not self.samtools_available and not self.warning_shown:
            print("\n" + "!"*70)
            print("WARNING: samtools not found in PATH")
            print("!"*70)
            print("BAM coverage and allele frequency will NOT be calculated.")
            print("VCF_AF will be used if available.")
            print("\nTo enable BAM processing:")
            print("  • Install samtools: https://www.htslib.org/download/")
            print("  • Add samtools to your PATH")
            print("  • On Windows: Add samtools.exe location to system PATH")
            print("!"*70 + "\n")
            self.warning_shown = True
    
    def get_coverage_and_af(self, bam_file: Path, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        if not self.samtools_available:
            return {'coverage': None, 'af': None, 'ref_count': None, 'alt_count': None}
        
        try:
            cmd = ['samtools', 'mpileup', '-r', f'{chrom}:{pos}-{pos}',
                   '-Q', '20', '-q', '20', str(bam_file)]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            
            if result.returncode == 0 and result.stdout.strip():
                fields = result.stdout.strip().split('\t')
                if len(fields) >= 5:
                    coverage = int(fields[3])
                    bases = fields[4]
                    
                    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                    i = 0
                    while i < len(bases):
                        if bases[i] == '^':
                            i += 2
                        elif bases[i] == '$':
                            i += 1
                        elif bases[i] in '.,':
                            base_counts[ref.upper()] += 1
                            i += 1
                        elif bases[i].upper() in 'ACGT':
                            base_counts[bases[i].upper()] += 1
                            i += 1
                        elif bases[i] in '+-':
                            i += 1
                            num = ''
                            while i < len(bases) and bases[i].isdigit():
                                num += bases[i]
                                i += 1
                            if num:
                                i += int(num)
                        else:
                            i += 1
                    
                    ref_count = base_counts.get(ref.upper(), 0)
                    alt_count = base_counts.get(alt.upper(), 0)
                    total = sum(base_counts.values())
                    af = alt_count / total if total > 0 else 0
                    
                    return {'coverage': coverage, 'af': af, 'ref_count': ref_count, 'alt_count': alt_count}
        except:
            pass
        
        return {'coverage': None, 'af': None, 'ref_count': None, 'alt_count': None}

# ============================================================================
# FILE MANAGEMENT
# ============================================================================

class Sample:
    """Container for sample files"""
    def __init__(self, name: str, vcf: Path, bam: Optional[Path] = None):
        self.name = name
        self.vcf = vcf
        self.bam = bam

class FileManagerCLI:
    """CLI mode: scan directory or single file"""
    
    def __init__(self, input_path: Path, gff_file: Path):
        self.input_path = input_path
        self.gff_file = gff_file
        self.samples = []
    
    def scan(self):
        print("\n" + "="*70)
        print("CLI MODE: SCANNING INPUT")
        print("="*70)
        
        if self.input_path.is_file():
            # Single VCF file
            if self.input_path.suffix == '.vcf':
                sample_name = self.input_path.stem
                self.samples.append(Sample(sample_name, self.input_path))
                print(f"Single VCF: {self.input_path.name}")
            else:
                raise ValueError(f"Input file must be .vcf: {self.input_path}")
        
        elif self.input_path.is_dir():
            # Directory of VCF files
            vcf_files = sorted(self.input_path.glob("*.vcf"))
            for vcf in vcf_files:
                sample_name = vcf.stem
                # Try to find matching BAM
                bam = self.input_path / f"{sample_name}.bam"
                self.samples.append(Sample(sample_name, vcf, bam if bam.exists() else None))
            
            print(f"VCF files found: {len(self.samples)}")
            for s in self.samples:
                bam_status = "✓ BAM" if s.bam else "✗ No BAM"
                print(f"  • {s.vcf.name} ({bam_status})")
        else:
            raise ValueError(f"Input path not found: {self.input_path}")
        
        print(f"GFF file: {self.gff_file.name}")
        print("="*70)
    
    def validate(self):
        if not self.samples:
            raise ValueError("No VCF files found")
        if not self.gff_file.exists():
            raise ValueError(f"GFF file not found: {self.gff_file}")

class FileManagerConfig:
    """Config mode: explicit sample-BAM pairing"""
    
    def __init__(self, config: Dict):
        self.config = config
        self.samples = []
        self.gff_file = None
    
    def scan(self):
        print("\n" + "="*70)
        print("CONFIG MODE: LOADING SAMPLES")
        print("="*70)
        
        # Get GFF
        if 'gff_file' not in self.config:
            raise ValueError("Config must specify 'gff_file'")
        self.gff_file = Path(self.config['gff_file'])
        
        # Get samples
        if 'samples' not in self.config:
            raise ValueError("Config must specify 'samples' list")
        
        for sample_config in self.config['samples']:
            name = sample_config.get('name')
            vcf_path = Path(sample_config['vcf'])
            bam_path = Path(sample_config['bam']) if 'bam' in sample_config else None
            
            if not name:
                name = vcf_path.stem
            
            self.samples.append(Sample(name, vcf_path, bam_path))
        
        print(f"Samples: {len(self.samples)}")
        for s in self.samples:
            bam_status = f"✓ BAM: {s.bam.name}" if s.bam else "✗ No BAM"
            print(f"  • {s.name}")
            print(f"    VCF: {s.vcf.name}")
            print(f"    {bam_status}")
        
        print(f"\nGFF file: {self.gff_file.name}")
        print("="*70)
    
    def validate(self):
        if not self.samples:
            raise ValueError("No samples in config")
        if not self.gff_file or not self.gff_file.exists():
            raise ValueError(f"GFF file not found: {self.gff_file}")
        
        for s in self.samples:
            if not s.vcf.exists():
                raise ValueError(f"VCF not found: {s.vcf}")
            if s.bam and not s.bam.exists():
                raise ValueError(f"BAM not found: {s.bam}")


# ============================================================================
# VARIANT PROCESSOR (keeping same as before)
# ============================================================================

class VariantProcessor:
    """Process variants with annotation"""
    
    def __init__(self, vcf_file: Path, gff_file: Path, bam_file: Optional[Path],
                 min_quality: float, min_coverage: int, color_mode: str, logger: Logger):
        self.vcf_file = vcf_file
        self.gff_file = gff_file
        self.bam_file = bam_file
        self.min_quality = min_quality
        self.min_coverage = min_coverage
        self.color_mode = color_mode
        self.logger = logger
        
        self.variants = []
        self.genes = []
        self.bam_parser = BAMParser()
        
        self.stats = {
            'total': 0,
            'passed': 0,
            'by_type': defaultdict(int),
            'by_location': defaultdict(int)
        }
    
    def load_genes(self):
        with open(self.gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                attrs = {}
                for attr in fields[8].split(';'):
                    if '=' in attr:
                        k, v = attr.split('=', 1)
                        attrs[k.strip()] = v.strip()
                
                name = attrs.get('Name') or attrs.get('gene_name') or attrs.get('gene') or attrs.get('ID') or 'Unknown'
                self.genes.append({
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'type': fields[2],
                    'strand': fields[6],
                    'name': name
                })
    
    def process(self):
        # Show samtools warning if BAM provided but samtools not available
        if self.bam_file and not self.bam_parser.samtools_available:
            self.bam_parser.show_warning()
        
        with open(self.vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                
                self.stats['total'] += 1
                
                chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
                qual = float(fields[5]) if fields[5] != '.' else None
                filter_status = fields[6] if len(fields) > 6 else '.'
                
                vcf_depth, vcf_af = None, None
                if len(fields) > 7:
                    info = fields[7]
                    dp_match = re.search(r'DP=(\d+)', info)
                    if dp_match:
                        vcf_depth = int(dp_match.group(1))
                    af_match = re.search(r'AF=([\d.]+)', info)
                    if af_match:
                        vcf_af = float(af_match.group(1))
                
                bam_data = {'coverage': None, 'af': None, 'ref_count': None, 'alt_count': None}
                if self.bam_file:
                    bam_data = self.bam_parser.get_coverage_and_af(self.bam_file, chrom, pos, ref, alt)
                    if self.stats['total'] % 100 == 0:
                        print(f"  Processed {self.stats['total']} variants...", end='\r')
                
                overlapping = [g for g in self.genes if g['chrom'] == chrom and g['start'] <= pos <= g['end']]
                location = self._get_location(overlapping)
                var_type = self._get_var_type(ref, alt)
                
                passes = True
                if self.min_quality > 0 and qual and qual < self.min_quality:
                    passes = False
                if self.min_coverage > 0 and bam_data['coverage'] and bam_data['coverage'] < self.min_coverage:
                    passes = False
                
                if passes:
                    self.stats['passed'] += 1
                
                self.stats['by_type'][var_type] += 1
                self.stats['by_location'][location] += 1
                
                self.variants.append({
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'qual': qual,
                    'filter': filter_status,
                    'vcf_depth': vcf_depth,
                    'vcf_af': vcf_af,
                    'bam_coverage': bam_data['coverage'],
                    'bam_af': bam_data['af'],
                    'bam_ref': bam_data['ref_count'],
                    'bam_alt': bam_data['alt_count'],
                    'genes': overlapping,
                    'location': location,
                    'type': var_type,
                    'passes': passes
                })
        
        print(f"  Processed {self.stats['total']} variants - Complete")
    
    def _get_location(self, genes: List[Dict]) -> str:
        if not genes:
            return 'intergenic'
        types = [g['type'] for g in genes]
        if 'CDS' in types:
            return 'coding'
        if 'exon' in types:
            return 'exon'
        if 'gene' in types:
            return 'gene'
        return 'genic'
    
    def _get_var_type(self, ref: str, alt: str) -> str:
        if len(ref) == 1 and len(alt) == 1:
            return 'SNV'
        elif len(ref) < len(alt):
            return 'INS'
        elif len(ref) > len(alt):
            return 'DEL'
        return 'MNV'
    
    def write_bed(self, bed_file: Path, sample_name: str):
        passed = [v for v in self.variants if v['passes']]
        
        with open(bed_file, 'w') as f:
            f.write(f'track name="{sample_name}" description="Annotated variants" itemRgb="On" visibility=2\n')
            
            for v in passed:
                start, end = v['pos'] - 1, v['pos']
                
                parts = [
                    f"Var:{v['ref']}>{v['alt']}",
                    f"Type:{v['type']}",
                ]
                
                if v['genes']:
                    parts.append(f"Gene:{','.join([g['name'] for g in v['genes'][:2]])}")
                else:
                    parts.append("Gene:intergenic")
                
                parts.append(f"Loc:{v['location']}")
                
                if v['qual']:
                    parts.append(f"QUAL:{v['qual']:.0f}")
                
                if v['bam_coverage']:
                    parts.append(f"Cov:{v['bam_coverage']}x")
                if v['bam_af'] is not None:
                    parts.append(f"BAM_AF:{v['bam_af']:.2f}")
                if v['vcf_af']:
                    parts.append(f"VCF_AF:{v['vcf_af']:.2f}")
                
                parts.append(f"Pos:{v['chrom']}:{v['pos']}")
                
                name = " | ".join(parts)
                rgb = self._get_color(v)
                score = self._get_score(v)
                strand = v['genes'][0]['strand'] if v['genes'] else '.'
                
                f.write(f"{v['chrom']}\t{start}\t{end}\t{name}\t{score}\t{strand}\t{start}\t{end}\t{rgb}\n")
    
    def _get_color(self, v: Dict) -> str:
        if self.color_mode == 'position':
            colors = {'coding': '255,0,0', 'exon': '255,128,0', 'gene': '0,128,255', 'genic': '0,255,128', 'intergenic': '128,128,128'}
            return colors.get(v['location'], '128,128,128')
        elif self.color_mode == 'type':
            colors = {'SNV': '0,0,255', 'INS': '0,255,0', 'DEL': '255,0,0', 'MNV': '255,0,255'}
            return colors.get(v['type'], '128,128,128')
        elif self.color_mode == 'coverage':
            cov = v['bam_coverage']
            if not cov:
                return '200,200,200'
            elif cov < 20:
                return '255,200,200'
            elif cov < 50:
                return '255,255,0'
            elif cov < 100:
                return '0,255,0'
            else:
                return '0,200,0'
        elif self.color_mode == 'allele_freq':
            af = v['bam_af'] if v['bam_af'] is not None else v['vcf_af']
            if af is None:
                return '200,200,200'
            elif af < 0.1:
                return '200,200,200'
            elif af < 0.25:
                return '150,150,255'
            elif af < 0.4:
                return '0,0,255'
            elif af < 0.6:
                return '0,150,0'
            else:
                return '255,0,0'
        return '128,128,128'
    
    def _get_score(self, v: Dict) -> int:
        base = {'coding': 1000, 'exon': 800, 'gene': 600, 'genic': 400, 'intergenic': 200}
        score = base.get(v['location'], 500)
        if v['qual']:
            score = int(score * min(v['qual'] / 100, 1.0))
        return score
    
    def write_tsv(self, tsv_file: Path):
        passed = [v for v in self.variants if v['passes']]
        
        with open(tsv_file, 'w') as f:
            headers = ['Chrom', 'Pos', 'Ref', 'Alt', 'Type', 'QUAL', 'Filter',
                      'Gene', 'Location', 'VCF_DP', 'VCF_AF', 'BAM_Cov', 'BAM_AF', 'BAM_Ref', 'BAM_Alt']
            f.write('\t'.join(headers) + '\n')
            
            for v in passed:
                gene = ','.join([g['name'] for g in v['genes']]) if v['genes'] else 'intergenic'
                row = [
                    v['chrom'], str(v['pos']), v['ref'], v['alt'], v['type'],
                    f"{v['qual']:.1f}" if v['qual'] else 'N/A',
                    v['filter'],
                    gene, v['location'],
                    str(v['vcf_depth']) if v['vcf_depth'] else 'N/A',
                    f"{v['vcf_af']:.3f}" if v['vcf_af'] else 'N/A',
                    str(v['bam_coverage']) if v['bam_coverage'] else 'N/A',
                    f"{v['bam_af']:.3f}" if v['bam_af'] is not None else 'N/A',
                    str(v['bam_ref']) if v['bam_ref'] is not None else 'N/A',
                    str(v['bam_alt']) if v['bam_alt'] is not None else 'N/A'
                ]
                f.write('\t'.join(row) + '\n')


# ============================================================================
# HEATMAP GENERATOR (with RED colors)
# ============================================================================

class HeatmapGenerator:
    """Generate heatmaps from TSV - RED color scheme - bedGraph only"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.heatmap_dir = output_dir / "heatmaps"
        self.heatmap_dir.mkdir(exist_ok=True)
    
    def generate_all(self, tsv_files: List[Path], gff_file: Path, window_size: int, with_bam: bool):
        print("\n" + "="*70)
        print("GENERATING HEATMAPS (RED COLOR SCHEME)")
        print("="*70)
        
        for tsv_file in tsv_files:
            sample = tsv_file.stem.replace('_details', '')
            print(f"\n[{sample}]")
            self._gene_hotspot(tsv_file, gff_file, sample)
            self._window_density(tsv_file, gff_file, sample, window_size)
            
            if with_bam:
                self._af_landscape(tsv_file, gff_file, sample, window_size)
        
        print("\n✓ Heatmap generation complete")
    
    def _gene_hotspot(self, tsv_file: Path, gff_file: Path, sample: str):
        print("  Gene hotspots...", end=' ')
        genes = self._load_genes(gff_file)
        variants = self._load_variants(tsv_file)
        
        gene_counts = {}
        for gene in genes:
            count = sum(1 for v in variants if v['chrom'] == gene['chrom'] and gene['start'] <= v['pos'] <= gene['end'])
            if count > 0:
                gene_counts[gene['name']] = {'chrom': gene['chrom'], 'start': gene['start'], 'end': gene['end'], 'count': count}
        
        bg_file = self.heatmap_dir / f"{sample}_gene_hotspots.bedGraph"
        self._write_bedgraph(gene_counts, bg_file, f"{sample}_genes", "Gene-level variant hotspots")
        print("✓")
    
    def _window_density(self, tsv_file: Path, gff_file: Path, sample: str, window: int):
        print(f"  Window density ({window}bp)...", end=' ')
        variants = self._load_variants(tsv_file)
        chrom_sizes = self._get_chrom_sizes(gff_file)
        
        windows = defaultdict(lambda: defaultdict(int))
        for v in variants:
            if v['chrom'] in chrom_sizes:
                w_start = (v['pos'] // window) * window
                windows[v['chrom']][w_start] += 1
        
        bg_file = self.heatmap_dir / f"{sample}_density_{window}bp.bedGraph"
        with open(bg_file, 'w') as f:
            # RED color scheme
            f.write(f'track type=bedGraph name="{sample}_density" ')
            f.write('description="Variant density" color=255,0,0 altColor=200,0,0\n')
            for chrom in sorted(windows):
                for w_start in sorted(windows[chrom]):
                    f.write(f"{chrom}\t{w_start}\t{w_start+window}\t{windows[chrom][w_start]}\n")
        
        print("✓")
    
    def _af_landscape(self, tsv_file: Path, gff_file: Path, sample: str, window: int):
        """Generate allele frequency landscape heatmap - RED color"""
        print(f"  AF landscape ({window}bp)...", end=' ')
        
        variants_with_af = []
        with open(tsv_file) as f:
            header = f.readline().strip().split('\t')
            
            # Try BAM_AF first, fallback to VCF_AF
            bam_af_idx = None
            vcf_af_idx = None
            try:
                bam_af_idx = header.index('BAM_AF')
            except ValueError:
                pass
            try:
                vcf_af_idx = header.index('VCF_AF')
            except ValueError:
                pass
            
            if bam_af_idx is None and vcf_af_idx is None:
                print("✗ (no AF data)")
                return
            
            # Determine which AF to use
            af_source = "BAM_AF" if bam_af_idx is not None else "VCF_AF"
            af_idx = bam_af_idx if bam_af_idx is not None else vcf_af_idx
            
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) > af_idx:
                    try:
                        af_str = fields[af_idx]
                        if af_str != 'N/A':
                            af = float(af_str)
                            if 0 <= af <= 1:
                                variants_with_af.append({
                                    'chrom': fields[0],
                                    'pos': int(fields[1]),
                                    'af': af
                                })
                    except ValueError:
                        continue
        
        if not variants_with_af:
            print(f"✗ (no valid {af_source} data)")
            return
        
        chrom_sizes = self._get_chrom_sizes(gff_file)
        
        window_af_sum = defaultdict(lambda: defaultdict(float))
        window_af_count = defaultdict(lambda: defaultdict(int))
        
        for v in variants_with_af:
            if v['chrom'] in chrom_sizes:
                w_start = (v['pos'] // window) * window
                window_af_sum[v['chrom']][w_start] += v['af']
                window_af_count[v['chrom']][w_start] += 1
        
        bg_file = self.heatmap_dir / f"{sample}_af_landscape.bedGraph"
        with open(bg_file, 'w') as f:
            # RED color scheme for AF
            f.write(f'track type=bedGraph name="{sample}_AF_landscape" ')
            f.write(f'description="Allele frequency landscape ({af_source})" ')
            f.write('color=255,0,0 altColor=200,0,0 viewLimits=0:100\n')
            
            for chrom in sorted(window_af_sum):
                for w_start in sorted(window_af_sum[chrom]):
                    count = window_af_count[chrom][w_start]
                    avg_af = window_af_sum[chrom][w_start] / count
                    scaled_af = avg_af * 100
                    f.write(f"{chrom}\t{w_start}\t{w_start+window}\t{scaled_af:.2f}\n")
        
        print(f"✓ (using {af_source})")
    
    def _load_genes(self, gff_file: Path) -> List[Dict]:
        genes = []
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 9 and fields[2] == 'gene':
                    attrs = {}
                    for a in fields[8].split(';'):
                        if '=' in a:
                            k, v = a.split('=', 1)
                            attrs[k.strip()] = v.strip()
                    name = attrs.get('Name') or attrs.get('ID') or 'Unknown'
                    genes.append({'chrom': fields[0], 'start': int(fields[3]), 'end': int(fields[4]), 'name': name})
        return genes
    
    def _load_variants(self, tsv_file: Path) -> List[Dict]:
        variants = []
        with open(tsv_file) as f:
            next(f)
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 2:
                    variants.append({'chrom': fields[0], 'pos': int(fields[1])})
        return variants
    
    def _get_chrom_sizes(self, gff_file: Path) -> Dict[str, int]:
        sizes = {}
        with open(gff_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    chrom, end = fields[0], int(fields[4])
                    if chrom not in sizes or end > sizes[chrom]:
                        sizes[chrom] = end
        return sizes
    
    def _write_bedgraph(self, data: Dict, bg_file: Path, track_name: str, description: str):
        with open(bg_file, 'w') as f:
            # RED color scheme
            f.write(f'track type=bedGraph name="{track_name}" ')
            f.write(f'description="{description}" color=255,0,0 altColor=200,0,0\n')
            for name, info in sorted(data.items(), key=lambda x: (x[1]['chrom'], x[1]['start'])):
                f.write(f"{info['chrom']}\t{info['start']}\t{info['end']}\t{info['count']}\n")

class SummaryWriter:
    @staticmethod
    def write_sample(sample: str, stats: Dict, filters: Dict, output_file: Path):
        with open(output_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write(f"SAMPLE: {sample}\n")
            f.write("="*70 + "\n\n")
            f.write("FILTERING:\n")
            f.write(f"  Min QUAL: {filters['min_quality']}\n")
            f.write(f"  Min Coverage: {filters['min_coverage']}\n\n")
            f.write("VARIANTS:\n")
            f.write(f"  Total: {stats['total']}\n")
            f.write(f"  Passed: {stats['passed']} ({stats['passed']/stats['total']*100:.1f}%)\n\n")
            f.write("BY TYPE:\n")
            for vtype, count in sorted(stats['by_type'].items(), key=lambda x: x[1], reverse=True):
                f.write(f"  {vtype}: {count}\n")
            f.write("\nBY LOCATION:\n")
            for loc, count in sorted(stats['by_location'].items(), key=lambda x: x[1], reverse=True):
                f.write(f"  {loc}: {count}\n")
    
    @staticmethod
    def write_combined(all_stats: Dict, filters: Dict, output_file: Path):
        with open(output_file, 'w') as f:
            f.write("="*70 + "\n")
            f.write("COMBINED ANALYSIS SUMMARY\n")
            f.write("="*70 + "\n\n")
            f.write(f"Samples: {len(all_stats)}\n\n")
            f.write("FILTERING:\n")
            f.write(f"  Min QUAL: {filters['min_quality']}\n")
            f.write(f"  Min Coverage: {filters['min_coverage']}\n\n")
            f.write("PER-SAMPLE STATS:\n")
            f.write("-"*70 + "\n")
            f.write(f"{'Sample':<20} {'Total':>10} {'Passed':>10} {'Coding':>10}\n")
            f.write("-"*70 + "\n")
            for sample, stats in sorted(all_stats.items()):
                coding = stats['by_location'].get('coding', 0)
                f.write(f"{sample:<20} {stats['total']:>10} {stats['passed']:>10} {coding:>10}\n")


# ============================================================================
# CONFIG FUNCTIONS
# ============================================================================

def load_config(config_file: Path) -> Dict:
    """Load configuration from JSON file"""
    with open(config_file) as f:
        config = json.load(f)
    return config

def save_config_template(output_file: Path):
    """Save a template configuration file"""
    template = {
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
    
    with open(output_file, 'w') as f:
        json.dump(template, f, indent=2)
    
    print("="*70)
    print("CONFIG TEMPLATE GENERATED")
    print("="*70)
    print(f"\nFile: {output_file}")
    print("\nConfig mode features:")
    print("  • Specify samples with explicit VCF-BAM pairing")
    print("  • BAM filename doesn't need to match VCF")
    print("  • Process multiple samples with different BAM names")
    print("\nEdit the config file and run with:")
    print(f"  python vcf_visualizer_v2.1.py --config {output_file}")
    print("="*70)

# ============================================================================
# MAIN
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="VCF Visualizer v2.1 - Complete Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
TWO MODES: CLI or CONFIG

==============================================================================
MODE 1: CLI - Quick Analysis
==============================================================================

Usage: python vcf_visualizer_v2.1.py -i INPUT -g GFF -o OUTPUT [OPTIONS]

INPUT can be:
  • Single VCF file: -i sample.vcf
  • Folder with VCFs: -i vcf_folder/
    (automatically finds matching .bam files in same folder)

Examples:

1. Single VCF file:
   python vcf_visualizer_v2.1.py -i data/sample.vcf -g genes.gff -o results

2. Folder with VCFs (auto-detect BAMs):
   python vcf_visualizer_v2.1.py -i vcf_folder/ -g genes.gff -o results
   
   vcf_folder/
   ├── sample1.vcf
   ├── sample1.bam  ← Auto-detected
   ├── sample2.vcf
   └── sample2.bam  ← Auto-detected

==============================================================================
MODE 2: CONFIG - Multi-Sample with Custom BAM Names
==============================================================================

Usage: python vcf_visualizer_v2.1.py --config config.json

Generate template:
  python vcf_visualizer_v2.1.py --generate-config

Config format:
{
  "gff_file": "annotations.gff",
  "output": "results",
  "samples": [
    {
      "name": "tumor1",
      "vcf": "vcfs/patient1_tumor.vcf",
      "bam": "bams/patient1_tumor_aligned.bam"  ← Different name OK!
    },
    {
      "name": "tumor2",
      "vcf": "vcfs/patient2_tumor.vcf",
      "bam": "bams/patient2_mapped_reads.bam"
    }
  ]
}

Benefits:
  • BAM names don't need to match VCF names
  • Explicit sample-to-BAM pairing
  • Process multiple samples in one run
  • Reproducible configuration

==============================================================================

OPTIONS:
  -q, --min-quality QUAL    Min QUAL score (default: 0)
  -c, --min-coverage COV    Min BAM coverage (default: 0)
  -m, --color-mode MODE     position|type|coverage|allele_freq
  -w, --window-size SIZE    Heatmap window size (default: 500)

HEATMAPS:
  All heatmaps use RED color scheme for better visibility
  • Gene hotspots: Variant count per gene (red gradient)
  • Window density: Variants per genomic window (red gradient)
  • AF landscape: Allele frequency distribution (red gradient, BAM required)
        """
    )
    
    # Mode selection
    mode_group = parser.add_mutually_exclusive_group()
    mode_group.add_argument('--config', type=Path, help='Config file (multi-sample mode)')
    mode_group.add_argument('--generate-config', action='store_true', 
                           help='Generate config template')
    
    # CLI mode arguments
    parser.add_argument('-i', '--input', type=Path, 
                       help='Input VCF file or folder (CLI mode)')
    parser.add_argument('-g', '--gff', type=Path, 
                       help='GFF annotation file (CLI mode, required with -i)')
    parser.add_argument('-o', '--output', type=Path, help='Output directory')
    
    # Processing options
    parser.add_argument('-q', '--min-quality', type=float, default=0, help='Min QUAL')
    parser.add_argument('-c', '--min-coverage', type=int, default=0, help='Min coverage')
    parser.add_argument('-m', '--color-mode', 
                       choices=['position', 'type', 'coverage', 'allele_freq'],
                       default='position', help='BED color mode')
    parser.add_argument('-w', '--window-size', type=int, default=500, help='Window size')
    parser.add_argument('--version', action='version', version=f'v{__version__}')
    
    args = parser.parse_args()
    
    # Handle config generation
    if args.generate_config:
        save_config_template(Path('vcf_visualizer_config.json'))
        return
    
    # Determine mode and load files
    if args.config:
        # CONFIG MODE
        if not args.config.exists():
            print(f"Error: Config file not found: {args.config}")
            sys.exit(1)
        
        config = load_config(args.config)
        file_manager = FileManagerConfig(config)
        
        # Get parameters from config, but CLI args take precedence
        output_dir = args.output if args.output else Path(config.get('output', 'results'))
        min_quality = args.min_quality if args.min_quality != 0 else config.get('min_quality', 0)
        min_coverage = args.min_coverage if args.min_coverage != 0 else config.get('min_coverage', 0)
        color_mode = args.color_mode if args.color_mode != 'position' else config.get('color_mode', 'position')
        window_size = args.window_size if args.window_size != 500 else config.get('window_size', 500)
        
        # Show overrides
        overrides = []
        if args.output:
            overrides.append(f"output: {args.output}")
        if args.min_quality != 0:
            overrides.append(f"min_quality: {args.min_quality}")
        if args.min_coverage != 0:
            overrides.append(f"min_coverage: {args.min_coverage}")
        if args.color_mode != 'position':
            overrides.append(f"color_mode: {args.color_mode}")
        if args.window_size != 500:
            overrides.append(f"window_size: {args.window_size}")
        
        if overrides:
            print("\n" + "="*70)
            print("CLI OVERRIDES APPLIED:")
            for override in overrides:
                print(f"  {override}")
            print("="*70)
        
    else:
        # CLI MODE
        if not args.input:
            parser.error("CLI mode requires -i INPUT and -g GFF")
        if not args.gff:
            parser.error("CLI mode requires -g GFF file")
        if not args.output:
            parser.error("Output directory (-o) is required")
        
        if not args.input.exists():
            parser.error(f"Input not found: {args.input}")
        if not args.gff.exists():
            parser.error(f"GFF file not found: {args.gff}")
        
        file_manager = FileManagerCLI(args.input, args.gff)
        output_dir = args.output
        min_quality = args.min_quality
        min_coverage = args.min_coverage
        color_mode = args.color_mode
        window_size = args.window_size
    
    # Setup output
    output_dir.mkdir(parents=True, exist_ok=True)
    bed_dir = output_dir / "bed_files"
    table_dir = output_dir / "tables"
    summary_dir = output_dir / "summaries"
    bed_dir.mkdir(exist_ok=True)
    table_dir.mkdir(exist_ok=True)
    summary_dir.mkdir(exist_ok=True)
    
    # Logger
    logger = Logger(output_dir)
    logger.log_param('mode', 'config' if args.config else 'cli')
    logger.log_param('min_quality', min_quality)
    logger.log_param('min_coverage', min_coverage)
    logger.log_param('color_mode', color_mode)
    logger.log_param('window_size', window_size)
    
    print("="*70)
    print(f"VCF VISUALIZER v{__version__}")
    print("="*70)
    
    # Scan files
    logger.log_step("File detection")
    file_manager.scan()
    file_manager.validate()
    
    # Validate filters AFTER scanning files
    has_bam = any(s.bam for s in file_manager.samples)
    if min_coverage > 0 and not has_bam:
        parser.error("Coverage filter (-c) requires BAM files")
    if color_mode in ['coverage', 'allele_freq'] and not has_bam:
        parser.error(f"Color mode '{color_mode}' requires BAM files")
    
    # Log files
    logger.log_file('gff', file_manager.gff_file)
    for sample in file_manager.samples:
        logger.log_file(f'vcf_{sample.name}', sample.vcf)
        if sample.bam:
            logger.log_file(f'bam_{sample.name}', sample.bam)
    
    # Process variants
    logger.log_step("Variant processing")
    print("\n" + "="*70)
    print("PROCESSING VARIANTS")
    print("="*70)
    
    all_stats = {}
    tsv_files = []
    has_bam = False
    
    for sample in file_manager.samples:
        print(f"\n[{sample.name}]")
        print(f"  VCF: {sample.vcf.name}")
        if sample.bam:
            print(f"  BAM: {sample.bam.name}")
            has_bam = True
        print(f"  Loading genes...")
        
        processor = VariantProcessor(
            sample.vcf, file_manager.gff_file, sample.bam,
            min_quality, min_coverage, color_mode, logger
        )
        
        processor.load_genes()
        print(f"  Processing variants...")
        processor.process()
        
        print(f"  Writing BED...")
        processor.write_bed(bed_dir / f"{sample.name}.bed", sample.name)
        
        print(f"  Writing TSV...")
        tsv_file = table_dir / f"{sample.name}_details.tsv"
        processor.write_tsv(tsv_file)
        tsv_files.append(tsv_file)
        
        print(f"  Writing summary...")
        SummaryWriter.write_sample(
            sample.name, processor.stats,
            {'min_quality': min_quality, 'min_coverage': min_coverage},
            summary_dir / f"{sample.name}_summary.txt"
        )
        
        all_stats[sample.name] = processor.stats
        print(f"  ✓ Complete ({processor.stats['passed']}/{processor.stats['total']} variants)")
    
    # Combined summary
    logger.log_step("Combined summary")
    SummaryWriter.write_combined(
        all_stats,
        {'min_quality': min_quality, 'min_coverage': min_coverage},
        output_dir / "combined_summary.txt"
    )
    
    # Generate heatmaps
    logger.log_step("Heatmap generation")
    heatmap_gen = HeatmapGenerator(output_dir)
    heatmap_gen.generate_all(tsv_files, file_manager.gff_file, window_size, has_bam)
    
    # Save log
    logger.save()
    
    # Final summary
    print("\n" + "="*70)
    print("✓ PIPELINE COMPLETE")
    print("="*70)
    print(f"\nOutput: {output_dir.absolute()}")
    print(f"\nProcessed {len(file_manager.samples)} sample(s):")
    for sample_name, stats in all_stats.items():
        print(f"  • {sample_name}: {stats['passed']}/{stats['total']} variants")
    
    print("\nGenerated:")
    print("  bed_files/   - IGV variants (color-coded)")
    print("  tables/      - TSV data for analysis")
    print("  summaries/   - Per-sample + combined statistics")
    if has_bam:
        print("  heatmaps/    - bedGraph files: Gene hotspots + Window density + AF landscape (RED)")
    else:
        print("  heatmaps/    - bedGraph files: Gene hotspots + Window density (RED)")
    print("  processing_log.json - Full reproducibility record")
    print("\nHeatmap format: bedGraph (RED gradient)")
    print("Load in IGV: File → Load from File → Select .bedGraph files")
    print("="*70)

if __name__ == '__main__':
    main()