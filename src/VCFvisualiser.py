#!/usr/bin/env python3
"""
VCF Visualizer - Comprehensive Variant Annotation and Visualization Tool
Version 1.0

A complete solution for generating annotated BED files from VCF, BAM, and GFF files
for visualization in IGV with multiple coloring modes and filtering options.
"""

import argparse
import sys
import json
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from datetime import datetime
import subprocess
import re
import hashlib


__version__ = "1.0.0"
__author__ = "Genomics Analysis Pipeline"


class ConfigLogger:
    """Log configuration and parameters for reproducibility"""
    
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.log_file = output_dir / "processing_log.json"
        self.config = {
            'version': __version__,
            'timestamp': datetime.now().isoformat(),
            'command_line': ' '.join(sys.argv),
            'parameters': {},
            'input_files': {},
            'processing_steps': [],
            'statistics': {}
        }
    
    def log_parameter(self, key: str, value):
        """Log a parameter"""
        self.config['parameters'][key] = str(value)
    
    def log_input_file(self, file_type: str, file_path: Path):
        """Log input file with checksum for verification"""
        if file_path.exists():
            checksum = self._calculate_checksum(file_path)
            self.config['input_files'][file_type] = {
                'path': str(file_path.absolute()),
                'name': file_path.name,
                'size_bytes': file_path.stat().st_size,
                'checksum_md5': checksum,
                'timestamp': datetime.fromtimestamp(file_path.stat().st_mtime).isoformat()
            }
    
    def log_step(self, step: str):
        """Log processing step"""
        self.config['processing_steps'].append({
            'step': step,
            'timestamp': datetime.now().isoformat()
        })
    
    def log_statistic(self, key: str, value):
        """Log a statistic"""
        # Convert sets to lists for JSON serialization
        if isinstance(value, dict):
            value = {k: list(v) if isinstance(v, set) else v for k, v in value.items()}
        self.config['statistics'][key] = value
    
    def _calculate_checksum(self, file_path: Path) -> str:
        """Calculate MD5 checksum of file"""
        md5 = hashlib.md5()
        with open(file_path, 'rb') as f:
            for chunk in iter(lambda: f.read(4096), b""):
                md5.update(chunk)
        return md5.hexdigest()
    
    def save(self):
        """Save log to file"""
        with open(self.log_file, 'w') as f:
            json.dump(self.config, f, indent=2)


class ContigChecker:
    """Check and validate contig names across files"""
    
    @staticmethod
    def extract_vcf_contigs(vcf_file: Path) -> set:
        """Extract contig names from VCF"""
        contigs = set()
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields:
                    contigs.add(fields[0])
        return contigs
    
    @staticmethod
    def extract_gff_contigs(gff_file: Path) -> set:
        """Extract contig names from GFF"""
        contigs = set()
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if fields:
                    contigs.add(fields[0])
        return contigs
    
    @staticmethod
    def check_compatibility(vcf_contigs: set, gff_contigs: set) -> Dict:
        """Check contig compatibility between VCF and GFF"""
        common = vcf_contigs & gff_contigs
        vcf_only = vcf_contigs - gff_contigs
        gff_only = gff_contigs - vcf_contigs
        
        return {
            'compatible': len(common) > 0,
            'common_contigs': common,
            'vcf_only': vcf_only,
            'gff_only': gff_only,
            'overlap_percentage': len(common) / len(vcf_contigs) * 100 if vcf_contigs else 0
        }


class FileDetector:
    """Detect and validate input files"""
    
    def __init__(self, input_dir: Path):
        self.input_dir = input_dir
        self.vcf_files = []
        self.bam_files = []
        self.gff_files = []
        
    def scan(self):
        """Scan directory for files"""
        if not self.input_dir.exists():
            raise FileNotFoundError(f"Input directory not found: {self.input_dir}")
        
        for file_path in self.input_dir.iterdir():
            if not file_path.is_file():
                continue
            
            name_lower = file_path.name.lower()
            
            if name_lower.endswith(('.vcf', '.vcf.gz')):
                self.vcf_files.append(file_path)
            elif name_lower.endswith('.bam'):
                self.bam_files.append(file_path)
            elif name_lower.endswith(('.gff', '.gff3', '.gtf')):
                self.gff_files.append(file_path)
        
        self.vcf_files.sort()
        self.bam_files.sort()
        self.gff_files.sort()
    
    def validate(self):
        """Validate required files present"""
        if not self.vcf_files:
            raise ValueError("No VCF files found in input directory")
        if not self.gff_files:
            raise ValueError("No GFF file found in input directory")
        if len(self.gff_files) > 1:
            print(f"Warning: Multiple GFF files found, using: {self.gff_files[0].name}")
    
    def print_summary(self):
        """Print file detection summary"""
        print("\n" + "=" * 70)
        print("FILE DETECTION")
        print("=" * 70)
        print(f"Input directory: {self.input_dir.absolute()}\n")
        print(f"VCF files found: {len(self.vcf_files)}")
        for vcf in self.vcf_files:
            print(f"  • {vcf.name}")
        print(f"\nGFF files found: {len(self.gff_files)}")
        for gff in self.gff_files:
            print(f"  • {gff.name}")
        print(f"\nBAM files found: {len(self.bam_files)}")
        for bam in self.bam_files:
            print(f"  • {bam.name}")
        print("=" * 70 + "\n")


class VariantProcessor:
    """Process variants with comprehensive annotation"""
    
    def __init__(self, vcf_file: Path, gff_file: Path, bam_file: Optional[Path] = None,
                 min_quality: float = 0, min_coverage: int = 0, 
                 color_mode: str = 'position', logger: Optional[ConfigLogger] = None):
        self.vcf_file = vcf_file
        self.gff_file = gff_file
        self.bam_file = bam_file
        self.min_quality = min_quality
        self.min_coverage = min_coverage
        self.color_mode = color_mode
        self.logger = logger
        
        self.variants = []
        self.genes = []
        self.samtools_available = self._check_samtools()
        
        self.stats = {
            'total_variants': 0,
            'passed_filters': 0,
            'failed_quality': 0,
            'failed_coverage': 0,
            'by_type': {},
            'by_location': {},
            'by_contig': {}
        }
    
    def _check_samtools(self) -> bool:
        """Check if samtools is available"""
        try:
            subprocess.run(['samtools', '--version'], capture_output=True, timeout=2)
            return True
        except:
            return False
    
    def load_genes(self):
        """Load gene annotations from GFF"""
        if self.logger:
            self.logger.log_step(f"Loading genes from {self.gff_file.name}")
        
        print(f"→ Loading genes from {self.gff_file.name}...")
        
        with open(self.gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                attributes = {}
                for attr in fields[8].split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attributes[key.strip()] = value.strip()
                
                gene_name = (attributes.get('Name') or 
                           attributes.get('gene_name') or 
                           attributes.get('gene') or
                           attributes.get('ID') or
                           attributes.get('locus_tag') or
                           'Unknown')
                
                self.genes.append({
                    'chrom': fields[0],
                    'start': int(fields[3]),
                    'end': int(fields[4]),
                    'type': fields[2],
                    'strand': fields[6],
                    'name': gene_name
                })
        
        print(f"  Loaded {len(self.genes)} gene features")
    
    def get_bam_data(self, chrom: str, pos: int, ref: str, alt: str) -> Dict:
        """Get comprehensive BAM data"""
        result = {
            'coverage': None,
            'ref_count': None,
            'alt_count': None,
            'bam_af': None,
            'base_quality': None,
            'mapping_quality': None,
            'available': False
        }
        
        if not self.bam_file or not self.samtools_available:
            return result
        
        try:
            cmd = ['samtools', 'mpileup', '-r', f'{chrom}:{pos}-{pos}',
                   '-Q', '20', '-q', '20', str(self.bam_file)]
            mpileup = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
            
            if mpileup.returncode == 0 and mpileup.stdout.strip():
                fields = mpileup.stdout.strip().split('\t')
                if len(fields) >= 5:
                    coverage = int(fields[3])
                    base_string = fields[4]
                    
                    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
                    
                    i = 0
                    while i < len(base_string):
                        char = base_string[i]
                        if char == '^':
                            i += 2
                            continue
                        elif char == '$':
                            i += 1
                            continue
                        elif char in '.,':
                            base_counts[ref.upper()] += 1
                        elif char.upper() in 'ACGT':
                            base_counts[char.upper()] += 1
                        elif char in '+-':
                            i += 1
                            if i < len(base_string):
                                num_str = ''
                                while i < len(base_string) and base_string[i].isdigit():
                                    num_str += base_string[i]
                                    i += 1
                                if num_str:
                                    i += int(num_str) - 1
                        i += 1
                    
                    ref_count = base_counts.get(ref.upper(), 0)
                    alt_count = base_counts.get(alt.upper(), 0)
                    total = sum(base_counts.values())
                    
                    bam_af = alt_count / total if total > 0 else None
                    
                    result = {
                        'coverage': coverage,
                        'ref_count': ref_count,
                        'alt_count': alt_count,
                        'bam_af': bam_af,
                        'base_quality': 20,  # Min quality used
                        'mapping_quality': 20,  # Min quality used
                        'available': True
                    }
        except:
            pass
        
        return result
    
    def parse_vcf_info(self, fields: List[str]) -> Dict:
        """Parse VCF fields"""
        info = {
            'chrom': fields[0],
            'pos': int(fields[1]),
            'id': fields[2] if fields[2] != '.' else None,
            'ref': fields[3],
            'alt': fields[4],
            'qual': None,
            'filter': fields[6] if len(fields) > 6 else '.',
            'vcf_depth': None,
            'vcf_af': None,
            'variant_type': None
        }
        
        # Quality
        if len(fields) > 5 and fields[5] != '.':
            try:
                info['qual'] = float(fields[5])
            except:
                pass
        
        # INFO field
        if len(fields) > 7:
            info_str = fields[7]
            dp_match = re.search(r'DP=(\d+)', info_str)
            if dp_match:
                info['vcf_depth'] = int(dp_match.group(1))
            
            af_match = re.search(r'AF=([\d.]+)', info_str)
            if af_match:
                info['vcf_af'] = float(af_match.group(1))
        
        # Variant type
        info['variant_type'] = self._classify_variant_type(info['ref'], info['alt'])
        
        return info
    
    def _classify_variant_type(self, ref: str, alt: str) -> str:
        """Classify variant type"""
        if len(ref) == 1 and len(alt) == 1:
            return 'SNV'
        elif len(ref) < len(alt):
            return 'INS'
        elif len(ref) > len(alt):
            return 'DEL'
        else:
            return 'MNV'
    
    def process_variants(self):
        """Process all variants"""
        if self.logger:
            self.logger.log_step(f"Processing variants from {self.vcf_file.name}")
        
        print(f"→ Processing {self.vcf_file.name}...")
        
        with open(self.vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 5:
                    continue
                
                self.stats['total_variants'] += 1
                
                # Parse VCF info
                var_info = self.parse_vcf_info(fields)
                
                # Get BAM data if available
                bam_data = None
                if self.bam_file:
                    bam_data = self.get_bam_data(
                        var_info['chrom'], var_info['pos'],
                        var_info['ref'], var_info['alt']
                    )
                    if self.stats['total_variants'] % 100 == 0:
                        print(f"  Processed {self.stats['total_variants']} variants...")
                
                # Find overlapping genes
                genes = self._find_overlapping_genes(var_info['chrom'], var_info['pos'])
                location_type = self._get_location_type(genes)
                
                # Apply filters
                passes = self._apply_filters(var_info, bam_data)
                
                if passes:
                    self.stats['passed_filters'] += 1
                else:
                    if var_info['qual'] and var_info['qual'] < self.min_quality:
                        self.stats['failed_quality'] += 1
                    if bam_data and bam_data['coverage'] and bam_data['coverage'] < self.min_coverage:
                        self.stats['failed_coverage'] += 1
                
                # Update statistics
                var_type = var_info['variant_type']
                self.stats['by_type'][var_type] = self.stats['by_type'].get(var_type, 0) + 1
                self.stats['by_location'][location_type] = self.stats['by_location'].get(location_type, 0) + 1
                self.stats['by_contig'][var_info['chrom']] = self.stats['by_contig'].get(var_info['chrom'], 0) + 1
                
                # Store variant
                variant = {
                    **var_info,
                    'genes': genes,
                    'location_type': location_type,
                    'bam_data': bam_data,
                    'passes_filter': passes
                }
                
                self.variants.append(variant)
        
        print(f"  Total: {self.stats['total_variants']}, Passed filters: {self.stats['passed_filters']}")
    
    def _find_overlapping_genes(self, chrom: str, pos: int) -> List[Dict]:
        """Find genes overlapping position"""
        overlapping = []
        for gene in self.genes:
            if gene['chrom'] == chrom and gene['start'] <= pos <= gene['end']:
                overlapping.append(gene)
        return overlapping
    
    def _get_location_type(self, genes: List[Dict]) -> str:
        """Get location type from genes"""
        if not genes:
            return 'intergenic'
        types = [g['type'] for g in genes]
        if 'CDS' in types:
            return 'coding'
        elif 'exon' in types:
            return 'exon'
        elif 'gene' in types:
            return 'gene'
        return 'genic'
    
    def _apply_filters(self, var_info: Dict, bam_data: Optional[Dict]) -> bool:
        """Apply quality and coverage filters"""
        # Quality filter
        if self.min_quality > 0:
            if not var_info['qual'] or var_info['qual'] < self.min_quality:
                return False
        
        # Coverage filter
        if self.min_coverage > 0 and bam_data:
            if not bam_data['coverage'] or bam_data['coverage'] < self.min_coverage:
                return False
        
        return True
    
    def get_statistics(self) -> Dict:
        """Get processing statistics"""
        return self.stats


class BEDGenerator:
    """Generate BED files with various coloring modes"""
    
    COLOR_SCHEMES = {
        # Position-based (default)
        'position': {
            'coding': '255,0,0',      # Red
            'exon': '255,128,0',      # Orange
            'gene': '0,128,255',      # Blue
            'genic': '0,255,128',     # Green
            'intergenic': '128,128,128'  # Gray
        },
        # Variant type based
        'type': {
            'SNV': '0,0,255',         # Blue
            'INS': '0,255,0',         # Green
            'DEL': '255,0,0',         # Red
            'MNV': '255,0,255'        # Magenta
        },
        # Allele frequency based (gradients)
        'allele_freq': {
            'very_low': '200,200,200',    # 0-0.1
            'low': '150,150,255',         # 0.1-0.25
            'medium': '0,0,255',          # 0.25-0.4
            'high': '0,150,0',            # 0.4-0.6
            'very_high': '255,0,0'        # 0.6-1.0
        },
        # Quality based
        'quality': {
            'low': '255,200,200',         # <20
            'medium': '255,255,0',        # 20-40
            'high': '0,255,0',            # 40-60
            'very_high': '0,200,0'        # >60
        }
    }
    
    @staticmethod
    def generate_bed(variants: List[Dict], output_file: Path, color_mode: str = 'position',
                    sample_name: str = "Sample"):
        """Generate BED file with specified coloring mode"""
        
        filtered_variants = [v for v in variants if v['passes_filter']]
        
        with open(output_file, 'w') as f:
            # Track header
            f.write(f'track name="{sample_name}" ')
            f.write(f'description="Variants colored by {color_mode}" ')
            f.write('visibility=2 itemRgb="On" useScore=1\n')
            
            for variant in filtered_variants:
                chrom = variant['chrom']
                start = variant['pos'] - 1  # BED is 0-based
                end = variant['pos']
                
                # Build comprehensive name
                name_parts = BEDGenerator._build_name_parts(variant, color_mode)
                name = " | ".join(name_parts)
                
                # Get color based on mode
                rgb = BEDGenerator._get_color(variant, color_mode)
                
                # Calculate score
                score = BEDGenerator._calculate_score(variant, color_mode)
                
                # Strand
                strand = variant['genes'][0]['strand'] if variant['genes'] else '.'
                
                # Write BED line
                f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\t")
                f.write(f"{start}\t{end}\t{rgb}\n")
    
    @staticmethod
    def _build_name_parts(variant: Dict, color_mode: str) -> List[str]:
        """Build name field parts"""
        parts = []
        
        # Variant
        parts.append(f"Var:{variant['ref']}>{variant['alt']}")
        
        # Type
        parts.append(f"Type:{variant['variant_type']}")
        
        # Gene
        if variant['genes']:
            gene_names = [g['name'] for g in variant['genes']]
            parts.append(f"Gene:{','.join(gene_names[:2])}")  # Max 2 genes
        else:
            parts.append("Gene:intergenic")
        
        # Location
        parts.append(f"Loc:{variant['location_type']}")
        
        # Quality
        if variant['qual']:
            parts.append(f"QUAL:{variant['qual']:.1f}")
        
        # Coverage and AF (if BAM available)
        if variant['bam_data'] and variant['bam_data']['available']:
            parts.append(f"Cov:{variant['bam_data']['coverage']}x")
            if variant['bam_data']['bam_af'] is not None:
                parts.append(f"BAM_AF:{variant['bam_data']['bam_af']:.2f}")
        
        # VCF AF
        if variant['vcf_af']:
            parts.append(f"VCF_AF:{variant['vcf_af']:.2f}")
        
        # Position
        parts.append(f"Pos:{variant['chrom']}:{variant['pos']}")
        
        return parts
    
    @staticmethod
    def _get_color(variant: Dict, color_mode: str) -> str:
        """Get RGB color based on mode"""
        if color_mode == 'position':
            return BEDGenerator.COLOR_SCHEMES['position'].get(
                variant['location_type'], '128,128,128'
            )
        
        elif color_mode == 'type':
            return BEDGenerator.COLOR_SCHEMES['type'].get(
                variant['variant_type'], '128,128,128'
            )
        
        elif color_mode == 'allele_freq':
            af = None
            if variant['bam_data'] and variant['bam_data']['bam_af'] is not None:
                af = variant['bam_data']['bam_af']
            elif variant['vcf_af']:
                af = variant['vcf_af']
            
            if af is None:
                return '128,128,128'
            elif af < 0.1:
                return BEDGenerator.COLOR_SCHEMES['allele_freq']['very_low']
            elif af < 0.25:
                return BEDGenerator.COLOR_SCHEMES['allele_freq']['low']
            elif af < 0.4:
                return BEDGenerator.COLOR_SCHEMES['allele_freq']['medium']
            elif af < 0.6:
                return BEDGenerator.COLOR_SCHEMES['allele_freq']['high']
            else:
                return BEDGenerator.COLOR_SCHEMES['allele_freq']['very_high']
        
        elif color_mode == 'quality':
            qual = variant['qual']
            if qual is None:
                return '128,128,128'
            elif qual < 20:
                return BEDGenerator.COLOR_SCHEMES['quality']['low']
            elif qual < 40:
                return BEDGenerator.COLOR_SCHEMES['quality']['medium']
            elif qual < 60:
                return BEDGenerator.COLOR_SCHEMES['quality']['high']
            else:
                return BEDGenerator.COLOR_SCHEMES['quality']['very_high']
        
        return '128,128,128'
    
    @staticmethod
    def _calculate_score(variant: Dict, color_mode: str) -> int:
        """Calculate BED score"""
        base_scores = {
            'coding': 1000,
            'exon': 800,
            'gene': 600,
            'genic': 400,
            'intergenic': 200
        }
        
        base = base_scores.get(variant['location_type'], 500)
        
        # Adjust by quality
        if variant['qual']:
            qual_factor = min(variant['qual'] / 100, 1.0)
        else:
            qual_factor = 0.5
        
        return int(base * qual_factor)


class SummaryWriter:
    """Write summary files"""
    
    @staticmethod
    def write_txt_summary(stats: Dict, vcf_file: Path, output_file: Path,
                         color_mode: str, filters: Dict):
        """Write text summary"""
        with open(output_file, 'w') as f:
            f.write("=" * 70 + "\n")
            f.write("VARIANT ANALYSIS SUMMARY\n")
            f.write("=" * 70 + "\n\n")
            
            f.write(f"Sample: {vcf_file.stem}\n")
            f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("FILTERING CRITERIA:\n")
            f.write("-" * 50 + "\n")
            f.write(f"Minimum QUAL: {filters['min_quality']}\n")
            f.write(f"Minimum Coverage: {filters['min_coverage']}\n")
            f.write(f"Color mode: {color_mode}\n\n")
            
            f.write("VARIANT COUNTS:\n")
            f.write("-" * 50 + "\n")
            f.write(f"Total variants: {stats['total_variants']}\n")
            f.write(f"Passed filters: {stats['passed_filters']} ")
            f.write(f"({stats['passed_filters']/stats['total_variants']*100:.1f}%)\n")
            f.write(f"Failed quality: {stats['failed_quality']}\n")
            f.write(f"Failed coverage: {stats['failed_coverage']}\n\n")
            
            f.write("BY VARIANT TYPE:\n")
            f.write("-" * 50 + "\n")
            for vtype, count in sorted(stats['by_type'].items(), key=lambda x: x[1], reverse=True):
                f.write(f"{vtype:10s}: {count:5d} ({count/stats['total_variants']*100:5.1f}%)\n")
            f.write("\n")
            
            f.write("BY LOCATION:\n")
            f.write("-" * 50 + "\n")
            for loc, count in sorted(stats['by_location'].items(), key=lambda x: x[1], reverse=True):
                f.write(f"{loc:15s}: {count:5d} ({count/stats['total_variants']*100:5.1f}%)\n")
            f.write("\n")
            
            f.write("BY CONTIG:\n")
            f.write("-" * 50 + "\n")
            for contig, count in sorted(stats['by_contig'].items(), key=lambda x: x[1], reverse=True):
                f.write(f"{contig:15s}: {count:5d} ({count/stats['total_variants']*100:5.1f}%)\n")
    
    @staticmethod
    def write_tsv_details(variants: List[Dict], output_file: Path):
        """Write detailed TSV file"""
        filtered = [v for v in variants if v['passes_filter']]
        
        with open(output_file, 'w') as f:
            # Header
            headers = [
                'Chrom', 'Pos', 'Ref', 'Alt', 'Type', 'QUAL', 'Filter',
                'Gene', 'Location', 'VCF_DP', 'VCF_AF',
                'BAM_Cov', 'BAM_AF', 'BAM_Ref', 'BAM_Alt'
            ]
            f.write('\t'.join(headers) + '\n')
            
            # Data
            for v in filtered:
                gene = ','.join([g['name'] for g in v['genes']]) if v['genes'] else 'intergenic'
                
                bam_cov = v['bam_data']['coverage'] if v['bam_data'] and v['bam_data']['available'] else 'N/A'
                bam_af = f"{v['bam_data']['bam_af']:.3f}" if v['bam_data'] and v['bam_data']['bam_af'] is not None else 'N/A'
                bam_ref = v['bam_data']['ref_count'] if v['bam_data'] and v['bam_data']['available'] else 'N/A'
                bam_alt = v['bam_data']['alt_count'] if v['bam_data'] and v['bam_data']['available'] else 'N/A'
                
                row = [
                    v['chrom'],
                    str(v['pos']),
                    v['ref'],
                    v['alt'],
                    v['variant_type'],
                    f"{v['qual']:.1f}" if v['qual'] else 'N/A',
                    v['filter'],
                    gene,
                    v['location_type'],
                    str(v['vcf_depth']) if v['vcf_depth'] else 'N/A',
                    f"{v['vcf_af']:.3f}" if v['vcf_af'] else 'N/A',
                    str(bam_cov),
                    bam_af,
                    str(bam_ref),
                    str(bam_alt)
                ]
                f.write('\t'.join(row) + '\n')


def main():
    parser = argparse.ArgumentParser(
        description="VCF Visualizer - Comprehensive Variant Annotation and Visualization",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
COMPREHENSIVE VARIANT VISUALIZATION TOOL

This tool generates annotated BED files from VCF, BAM, and GFF files for
visualization in IGV with multiple coloring modes and filtering options.

INPUT:
  -i: Folder containing:
      • 1 or multiple VCF files
      • 1 GFF file (contig names checked)
      • BAM files (optional, matched by name)

OUTPUT:
  bed_files/     - Annotated BED files for IGV
  summaries/     - Text summaries for each sample
  tables/        - Detailed TSV tables
  processing_log.json - Complete processing log

COLOR MODES (default: position):
  position      - Color by genomic position (coding/exon/gene/intergenic)
  type          - Color by variant type (SNV/INS/DEL/MNV)
  allele_freq   - Color by allele frequency (requires -b)
  quality       - Color by quality evidence (requires -b)

EXAMPLES:

1. Default processing (multiple VCFs, position coloring):
   python vcf_visualizer.py -i input_folder -o results

2. Single VCF with BAM, color by allele frequency:
   python vcf_visualizer.py -i data -o results -b sample.bam -m allele_freq

3. Quality filtering with coverage filter:
   python vcf_visualizer.py -i data -o results -q 20 -c 30

4. Custom reference name:
   python vcf_visualizer.py -i data -o results -r hg38

REPRODUCIBILITY OPTIONS:
  --version              Show version
  --save-params FILE     Save parameters to file
  --load-params FILE     Load parameters from file
  -r, --reference        Reference genome name (for documentation)
  --check-contigs        Verify contig compatibility
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input-dir', required=True, type=Path,
                       help='Input directory with VCF, GFF, and optionally BAM files')
    parser.add_argument('-o', '--output-dir', required=True, type=Path,
                       help='Output directory')
    
    # Optional arguments
    parser.add_argument('-b', '--bam', type=Path,
                       help='Process single VCF with this BAM file (enables advanced color modes)')
    parser.add_argument('-q', '--min-quality', type=float, default=0,
                       help='Minimum QUAL score (default: 0)')
    parser.add_argument('-c', '--min-coverage', type=int, default=0,
                       help='Minimum coverage (requires -b, default: 0)')
    parser.add_argument('-m', '--color-mode', 
                       choices=['position', 'type', 'allele_freq', 'quality'],
                       default='position',
                       help='Color mode (default: position)')
    parser.add_argument('-r', '--reference', default='unspecified',
                       help='Reference genome name (for documentation)')
    
    # Reproducibility options
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('--check-contigs', action='store_true',
                       help='Check contig name compatibility between VCF and GFF')
    parser.add_argument('--save-params', type=Path,
                       help='Save parameters to JSON file')
    parser.add_argument('--load-params', type=Path,
                       help='Load parameters from JSON file')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.color_mode in ['allele_freq', 'quality'] and not args.bam:
        parser.error(f"Color mode '{args.color_mode}' requires -b/--bam option")
    
    if args.min_coverage > 0 and not args.bam:
        parser.error("Coverage filter requires -b/--bam option")
    
    # Create output directories
    args.output_dir.mkdir(parents=True, exist_ok=True)
    bed_dir = args.output_dir / 'bed_files'
    summary_dir = args.output_dir / 'summaries'
    table_dir = args.output_dir / 'tables'
    
    bed_dir.mkdir(exist_ok=True)
    summary_dir.mkdir(exist_ok=True)
    table_dir.mkdir(exist_ok=True)
    
    # Initialize logger
    logger = ConfigLogger(args.output_dir)
    logger.log_parameter('input_dir', args.input_dir)
    logger.log_parameter('output_dir', args.output_dir)
    logger.log_parameter('min_quality', args.min_quality)
    logger.log_parameter('min_coverage', args.min_coverage)
    logger.log_parameter('color_mode', args.color_mode)
    logger.log_parameter('reference', args.reference)
    
    print("=" * 70)
    print("VCF VISUALIZER v" + __version__)
    print("=" * 70)
    print()
    
    # Detect files
    print("STEP 1: FILE DETECTION")
    detector = FileDetector(args.input_dir)
    detector.scan()
    detector.validate()
    detector.print_summary()
    
    # Log input files
    for vcf in detector.vcf_files:
        logger.log_input_file(f'vcf_{vcf.stem}', vcf)
    logger.log_input_file('gff', detector.gff_files[0])
    if args.bam:
        logger.log_input_file('bam', args.bam)
    
    # Check contigs if requested
    if args.check_contigs:
        print("STEP 2: CONTIG COMPATIBILITY CHECK")
        print("=" * 70)
        vcf_contigs = ContigChecker.extract_vcf_contigs(detector.vcf_files[0])
        gff_contigs = ContigChecker.extract_gff_contigs(detector.gff_files[0])
        compat = ContigChecker.check_compatibility(vcf_contigs, gff_contigs)
        
        print(f"VCF contigs: {len(vcf_contigs)}")
        print(f"GFF contigs: {len(gff_contigs)}")
        print(f"Common contigs: {len(compat['common_contigs'])}")
        print(f"Overlap: {compat['overlap_percentage']:.1f}%")
        
        if not compat['compatible']:
            print("\n⚠️  WARNING: No common contigs found between VCF and GFF!")
            print("VCF contigs:", sorted(list(vcf_contigs))[:5], "...")
            print("GFF contigs:", sorted(list(gff_contigs))[:5], "...")
        else:
            print("✓ Compatible\n")
        
        logger.log_statistic('contig_check', compat)
    
    # Process files
    print("\nSTEP 3: VARIANT PROCESSING")
    print("=" * 70)
    
    gff_file = detector.gff_files[0]
    
    # Determine processing mode
    if args.bam:
        # Single file mode with BAM
        if len(detector.vcf_files) > 1:
            print(f"⚠️  BAM specified but multiple VCFs found. Processing first VCF only.")
        
        vcf_file = detector.vcf_files[0]
        processor = VariantProcessor(
            vcf_file, gff_file, args.bam,
            args.min_quality, args.min_coverage,
            args.color_mode, logger
        )
        processor.load_genes()
        processor.process_variants()
        
        # Generate outputs
        sample_name = vcf_file.stem
        BEDGenerator.generate_bed(
            processor.variants,
            bed_dir / f"{sample_name}.bed",
            args.color_mode,
            sample_name
        )
        
        SummaryWriter.write_txt_summary(
            processor.get_statistics(),
            vcf_file,
            summary_dir / f"{sample_name}_summary.txt",
            args.color_mode,
            {'min_quality': args.min_quality, 'min_coverage': args.min_coverage}
        )
        
        SummaryWriter.write_tsv_details(
            processor.variants,
            table_dir / f"{sample_name}_details.tsv"
        )
        
        logger.log_statistic(f'{sample_name}_stats', processor.get_statistics())
    
    else:
        # Multiple file mode (default: position coloring)
        for vcf_file in detector.vcf_files:
            processor = VariantProcessor(
                vcf_file, gff_file, None,
                args.min_quality, 0,
                'position', logger
            )
            processor.load_genes()
            processor.process_variants()
            
            # Generate outputs
            sample_name = vcf_file.stem
            BEDGenerator.generate_bed(
                processor.variants,
                bed_dir / f"{sample_name}.bed",
                'position',
                sample_name
            )
            
            SummaryWriter.write_txt_summary(
                processor.get_statistics(),
                vcf_file,
                summary_dir / f"{sample_name}_summary.txt",
                'position',
                {'min_quality': args.min_quality, 'min_coverage': 0}
            )
            
            SummaryWriter.write_tsv_details(
                processor.variants,
                table_dir / f"{sample_name}_details.tsv"
            )
            
            logger.log_statistic(f'{sample_name}_stats', processor.get_statistics())
    
    # Save log
    logger.save()
    
    # Final summary
    print("\n" + "=" * 70)
    print("✓ PROCESSING COMPLETE")
    print("=" * 70)
    print(f"\nOutput directory: {args.output_dir.absolute()}")
    print(f"\nGenerated files:")
    print(f"  BED files: {bed_dir}/")
    print(f"  Summaries: {summary_dir}/")
    print(f"  Tables: {table_dir}/")
    print(f"  Log: processing_log.json")
    print("\nTo view in IGV:")
    print(f"  Load BED files from: {bed_dir}/")
    print("=" * 70)


if __name__ == '__main__':
    main()