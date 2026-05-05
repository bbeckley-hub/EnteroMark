#!/usr/bin/env python3
"""
EnteroMark AMRfinderPlus - Enterococcus faecium Antimicrobial Resistance
Comprehensive AMR analysis with interactive per‑genome reports and clean batch summary
Author: Brown Beckley
Affiliation: University of Ghana Medical School
Version: 1.0.0 (Dynamic Database)
Uses BUNDLED AMRFinderPlus 4.2.7 with LATEST dynamic database
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Optional
import argparse
import re
from datetime import datetime
import psutil
import json
import random
from collections import defaultdict

class EnteroMarkAMRfinder:
    def __init__(self, cpus: int = None):
        # Setup logging FIRST
        self.logger = self._setup_logging()
        
        # Get module directory and set bundled paths
        self.module_dir = os.path.dirname(os.path.abspath(__file__))
        self.bundled_amrfinder = os.path.join(self.module_dir, "bin", "amrfinder")
        self.bundled_update = os.path.join(self.module_dir, "bin", "amrfinder_update")
        
        # Initialize available_ram before calculating cpus
        self.available_ram = self._get_available_ram()
        
        # Then calculate resources - MAXIMUM SPEED MODE
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # DYNAMIC DATABASE: find the latest dated folder (starts with 20)
        self.bundled_database = self._get_latest_database()
        
        # If no database found, log warning (analysis will be skipped later)
        if self.bundled_database is None:
            self.logger.warning("No AMRfinderPlus database found. Please run with --update-db to download.")
        
        # Read database version or set to Unknown
        db_version = self._get_database_version() if self.bundled_database else "Unknown"
        
        # ==================== E. FAECIUM SPECIFIC GENE SETS ====================
        # CRITICAL Vancomycin resistance (high-level)
        self.critical_vancomycin = {
            'vanA', 'vanB', 'vanD', 'vanM', 'vanZ-A', 'vanY-A', 'vanX-A', 'vanA', 'vanH-A', 'vanS-A', 'vanR-A'
        }
        # CRITICAL Linezolid resistance
        self.critical_linezolid = {
            'optrA', 'cfr', 'cfrB', 'poxtA'
        }
        # CRITICAL Aminoglycoside resistance 
        self.critical_aminoglycoside = {
            'aac(6\')-Ii', 'ant(6)-Ia', 'aph(3\')-III', 'aac(6\')-Ie-aph(2\'\')-Ia'
        }
        
        # HIGH RISK resistance genes for E. faecium 
        self.high_risk_resistance = {
            # Vancomycin low-level / rare types
            'vanC', 'vanE', 'vanG', 'vanL', 'vanN','vanA', 'vanB', 'vanD', 'vanM', 'vanZ-A', 'vanY-A', 'vanX-A', 'vanA', 'vanH-A', 'vanS-A', 'vanR-A'
            # Macrolide-lincosamide-streptogramin B
            'ermA', 'ermB', 'ermC', 'ermF', 'ermX', 'msrA', 'msrB', 'mefA', 'mphA',
            # Tetracycline
            'tetM', 'tetL', 'tetO', 'tetK', 'tetS',
            # Aminoglycoside (non-critical)
            'aac(6\')-Ii', 'ant(6)-Ia', 'aph(3\')-III', 'aac(6\')-Ie-aph(2\'\')-Ia',   
            # Fluoroquinolone (mostly point mutations, but some genes)
            'qnrA', 'qnrB', 'qnrS',
            # Chloramphenicol
            'cat', 'catA', 'catB',
            # Beta-lactam 
            'pbp5',
            # Efflux pumps (overexpression contributes)
            'msrC', 'emeA', 'efmA',
            # Other
            'dfrG',  # trimethoprim
            'sul1', 'sul2'  # sulfonamide
        }
        
        # Union of all critical risk genes
        self.critical_risk_genes = set.union(
            self.critical_vancomycin,
            self.critical_linezolid,
            self.critical_aminoglycoside
        )
        # Union of all high-risk genes
        self.high_risk_genes = set.union(
            self.critical_risk_genes,
            self.high_risk_resistance
        )
        
        # ASCII art with warm amber/orange theme
        self.ascii_art = r"""
╔═══════════════════════════════════════════════════════════════════════════════╗
║                                                                               ║
║   ███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███╗   ███╗ █████╗ ██████╗ ██╗  ██╗
║   ██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗████╗ ████║██╔══██╗██╔══██╗██║ ██╔╝
║   █████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║██╔████╔██║███████║██████╔╝█████╔╝ 
║   ██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║██║╚██╔╝██║██╔══██║██╔══██╗██╔═██╗ 
║   ███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝██║ ╚═╝ ██║██║  ██║██║  ██║██║  ██╗
║   ╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝
║                                                                               ║
║                     EnteroMark - E. faecium AMR Analysis                      ║
║                                                                               ║
╚═══════════════════════════════════════════════════════════════════════════════╝
"""
        self.metadata = {
            "tool_name": "EnteroMark AMRfinderPlus",
            "version": "1.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "amrfinder_version": "4.2.7",
            "database_version": db_version
        }
        
        self.science_quotes = [
            "The important thing is not to stop questioning. Curiosity has its own reason for existing. - Albert Einstein",
            "Science is not only a disciple of reason but also one of romance and passion. - Stephen Hawking",
            "Somewhere, something incredible is waiting to be known. - Carl Sagan",
            "The good thing about science is that it's true whether or not you believe in it. - Neil deGrasse Tyson",
            "In science, there are no shortcuts to truth. - Karl Popper",
            "Science knows no country, because knowledge belongs to humanity. - Louis Pasteur",
            "The science of today is the technology of tomorrow. - Edward Teller",
            "Nothing in life is to be feared, it is only to be understood. - Marie Curie",
            "EnteroMark turns genomic complexity into actionable insights for VRE surveillance. - Brown Beckley"
        ]
    
    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> float:
        try:
            ram_gb = psutil.virtual_memory().available / (1024 ** 3)
            return ram_gb
        except Exception:
            return 8  # fallback
    
    def _calculate_optimal_cpus(self, user_cpus: int = None) -> int:
        if user_cpus is not None:
            self._log_resource_info(user_cpus)
            return user_cpus
        try:
            total_physical_cores = psutil.cpu_count(logical=False) or os.cpu_count() or 2
            if total_physical_cores <= 4:
                optimal_cpus = total_physical_cores
            elif total_physical_cores <= 8:
                optimal_cpus = total_physical_cores - 1
            elif total_physical_cores <= 16:
                optimal_cpus = max(8, total_physical_cores - 2)
            elif total_physical_cores <= 32:
                optimal_cpus = max(16, total_physical_cores - 3)
            else:
                optimal_cpus = min(32, int(total_physical_cores * 0.95))
            optimal_cpus = max(1, min(optimal_cpus, total_physical_cores))
            self._log_resource_info(optimal_cpus, total_physical_cores)
            return optimal_cpus
        except Exception:
            return os.cpu_count() or 4
    
    def _calculate_max_concurrent(self, genome_count: int) -> int:
        """RAM-aware concurrency: limit jobs based on memory per genome (1.5 GB assumed)"""
        mem_per_job = 1.5  # GB, typical for AMRfinder on E. faecium
        max_by_ram = max(1, int(self.available_ram / mem_per_job))
        max_by_cpu = self.cpus
        concurrent = min(max_by_cpu, max_by_ram, genome_count)
        self.logger.info(f"RAM-aware concurrency: {concurrent} jobs (RAM: {self.available_ram:.1f} GB, CPUs: {self.cpus})")
        return concurrent
    
    def _log_resource_info(self, cpus: int, total_cores: int = None):
        self.logger.info(f"Available RAM: {self.available_ram:.1f} GB")
        if total_cores:
            self.logger.info(f"System CPU cores: {total_cores}")
            utilization = (cpus / total_cores) * 100
            self.logger.info(f"Using CPU cores: {cpus} ({utilization:.1f}% of available cores)")
        else:
            self.logger.info(f"Using user-specified CPU cores: {cpus}")
        if cpus <= 4:
            self.logger.info("💡 Performance: Multi-core (max speed for small systems)")
        elif cpus <= 8:
            self.logger.info("💡 Performance: High-speed mode")
        else:
            self.logger.info("💡 Performance: MAXIMUM SPEED MODE 🚀")
    
    def _get_latest_database(self) -> Optional[str]:
        """Find the latest dated database folder in data/amrfinder_db/ (starts with 20)"""
        db_root = os.path.join(self.module_dir, "data", "amrfinder_db")
        if not os.path.exists(db_root):
            self.logger.warning(f"Database root directory not found: {db_root}")
            return None
        # Find all subdirectories starting with '20'
        candidates = []
        for item in os.listdir(db_root):
            full_path = os.path.join(db_root, item)
            if os.path.isdir(full_path) and item.startswith('20'):
                candidates.append(item)
        if not candidates:
            self.logger.warning("No database folder starting with '20' found.")
            return None
        # Sort lexicographically (YYYY-MM-DD works) and take the latest
        latest = sorted(candidates)[-1]
        latest_path = os.path.join(db_root, latest)
        self.logger.info(f"Using latest database: {latest_path}")
        return latest_path
    
    def _get_database_version(self) -> str:
        """Read version.txt from the database folder or fallback to folder name"""
        if not self.bundled_database:
            return "Unknown"
        version_file = os.path.join(self.bundled_database, "version.txt")
        if os.path.exists(version_file):
            with open(version_file, 'r') as f:
                return f.read().strip()
        # Fallback to folder name
        return os.path.basename(self.bundled_database)
    
    def update_database(self) -> bool:
        """Download the latest AMRfinderPlus database using bundled amrfinder_update"""
        if not os.path.exists(self.bundled_update):
            self.logger.error(f"amrfinder_update not found at {self.bundled_update}")
            return False
        if not os.access(self.bundled_update, os.X_OK):
            self.logger.warning("amrfinder_update not executable, fixing permissions...")
            os.chmod(self.bundled_update, 0o755)
        db_dir = os.path.join(self.module_dir, "data", "amrfinder_db")
        os.makedirs(db_dir, exist_ok=True)
        self.logger.info("Updating AMRfinderPlus database...")
        try:
            cmd = [self.bundled_update, "--database", db_dir]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info("Database update completed successfully.")
            # Re‑detect latest database
            self.bundled_database = self._get_latest_database()
            if self.bundled_database:
                self.metadata['database_version'] = self._get_database_version()
                self.logger.info(f"New database version: {self.metadata['database_version']}")
                return True
            else:
                self.logger.error("Database update succeeded but no database folder found.")
                return False
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Database update failed: {e}")
            self.logger.error(f"STDERR: {e.stderr}")
            return False
    
    def check_amrfinder_installed(self) -> bool:
        try:
            if not os.path.exists(self.bundled_amrfinder):
                self.logger.error(f"Bundled AMRfinderPlus not found at: {self.bundled_amrfinder}")
                return False
            if not os.access(self.bundled_amrfinder, os.X_OK):
                self.logger.warning("Fixing permissions on bundled AMRfinderPlus")
                os.chmod(self.bundled_amrfinder, 0o755)
            result = subprocess.run([self.bundled_amrfinder, '--version'],
                                    capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info(f"Bundled AMRfinderPlus version: {version_line}")
            
            # Check dynamic database
            if self.bundled_database and os.path.exists(self.bundled_database):
                self.logger.info(f"✅ Bundled database found: {self.bundled_database}")
                db_version_file = os.path.join(self.bundled_database, "version.txt")
                if os.path.exists(db_version_file):
                    with open(db_version_file, 'r') as f:
                        db_version = f.read().strip()
                        self.logger.info(f"✅ Database version: {db_version}")
                else:
                    self.logger.info(f"✅ Database folder: {os.path.basename(self.bundled_database)}")
                return True
            else:
                self.logger.warning(f"⚠️ Bundled database not found at expected location.")
                self.logger.info("Please run with --update-db to download the latest database.")
                return False
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(f"Bundled AMRfinderPlus check failed: {e}")
            return False
    
    def run_amrfinder_single_genome(self, genome_file: str, output_dir: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"{genome_name}_amrfinder.txt")
        run_cpus = self.cpus
        
        cmd = [
            self.bundled_amrfinder,
            '-n', genome_file,
            '--output', output_file,
            '--threads', str(run_cpus),
            '--plus',
            '--organism', 'Enterococcus_faecium'
        ]
        if self.bundled_database and os.path.exists(self.bundled_database):
            cmd.extend(['--database', self.bundled_database])
            self.logger.info(f"Using bundled database: {self.bundled_database}")
        else:
            self.logger.warning("No database found, using default AMRfinderPlus location (may fail).")
        
        self.logger.info(f"Running BUNDLED AMRfinderPlus: {genome_name} (using {run_cpus} CPU cores)")
        try:
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            hits = self._parse_amrfinder_output(output_file)
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
        except subprocess.CalledProcessError as e:
            self.logger.error(f"AMRfinderPlus failed for {genome_name}: {e.stderr}")
            return {
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_amrfinder_output(self, amrfinder_file: str) -> List[Dict]:
        hits = []
        try:
            with open(amrfinder_file, 'r') as f:
                lines = f.readlines()
            if not lines or len(lines) < 2:
                return hits
            headers = lines[0].strip().split('\t')
            for line_num, line in enumerate(lines[1:], 2):
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) >= len(headers):
                    hit = dict(zip(headers, parts))
                    # Keep all original fields, but also add standardised ones
                    processed_hit = {
                        'protein_id': hit.get('Protein id', ''),
                        'contig_id': hit.get('Contig id', ''),
                        'start': hit.get('Start', ''),
                        'stop': hit.get('Stop', ''),
                        'strand': hit.get('Strand', ''),
                        'gene_symbol': str(hit.get('Element symbol', '')),
                        'sequence_name': hit.get('Element name', ''),
                        'scope': hit.get('Scope', ''),
                        'element_type': hit.get('Type', ''),
                        'element_subtype': hit.get('Subtype', ''),
                        'class': hit.get('Class', ''),
                        'subclass': hit.get('Subclass', ''),
                        'method': hit.get('Method', ''),
                        'target_length': hit.get('Target length', ''),
                        'ref_length': hit.get('Reference sequence length', ''),
                        'coverage': hit.get('% Coverage of reference', '').replace('%', ''),
                        'identity': hit.get('% Identity to reference', '').replace('%', ''),
                        'alignment_length': hit.get('Alignment length', ''),
                        'accession': hit.get('Closest reference accession', ''),
                        'closest_name': hit.get('Closest reference name', ''),
                        'hmm_id': hit.get('HMM accession', ''),
                        'hmm_description': hit.get('HMM description', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning(f"Line {line_num} has {len(parts)} parts, expected {len(headers)}")
        except Exception as e:
            self.logger.error(f"Error parsing {amrfinder_file}: {e}")
        self.logger.info(f"Parsed {len(hits)} AMR hits from {amrfinder_file}")
        return hits
    
    def _create_amrfinder_html_report(self, genome_name: str, hits: List[Dict], output_dir: str):
        """Interactive per‑genome HTML report with ALL columns from AMRfinder, scrollable tables, warm theme."""
        analysis = self._analyze_efaecalis_amr_results(hits)
        
        # JavaScript for interactive features (search, export, print)
        interactive_js = f"""
        <script>
            function searchTable(tableId, searchTerm) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let visibleCount = 0;
                for (let i = 1; i < rows.length; i++) {{
                    const row = rows[i];
                    const text = row.textContent.toLowerCase();
                    if (text.includes(searchTerm.toLowerCase())) {{
                        row.style.display = '';
                        visibleCount++;
                    }} else {{
                        row.style.display = 'none';
                    }}
                }}
                const resultCounter = document.getElementById('result-counter-' + tableId);
                if (resultCounter) resultCounter.textContent = visibleCount + ' results found';
            }}
            function exportToCSV(tableId, filename) {{
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                let csv = [];
                const headerCells = rows[0].getElementsByTagName('th');
                const headerRow = [];
                for (let cell of headerCells) headerRow.push(cell.textContent);
                csv.push(headerRow.join(','));
                for (let i = 1; i < rows.length; i++) {{
                    if (rows[i].style.display !== 'none') {{
                        const cells = rows[i].getElementsByTagName('td');
                        const row = [];
                        for (let cell of cells) {{
                            let text = cell.textContent.trim();
                            text = text.replace(/\\n/g, ' ').replace(/,/g, ';');
                            row.push(text);
                        }}
                        csv.push(row.join(','));
                    }}
                }}
                const blob = new Blob([csv.join('\\n')], {{ type: 'text/csv' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function exportToJSON(dataVar, filename) {{
                const data = window[dataVar];
                const jsonStr = JSON.stringify(data, null, 2);
                const blob = new Blob([jsonStr], {{ type: 'application/json' }});
                const url = window.URL.createObjectURL(blob);
                const a = document.createElement('a');
                a.href = url;
                a.download = filename;
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
                window.URL.revokeObjectURL(url);
            }}
            function printReport() {{
                const printWindow = window.open('', '_blank');
                printWindow.document.write('<html><head><title>{genome_name} - AMR Report</title><style>body {{ font-family: Arial; margin: 20px; }} .no-print {{ display: none; }} table {{ border-collapse: collapse; width: 100%; }} th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }} .critical-row {{ background-color: #f8d7da; }} .high-risk-row {{ background-color: #fff3cd; }}</style></head><body>');
                const content = document.querySelector('.container').cloneNode(true);
                const noPrint = content.querySelectorAll('.no-print');
                noPrint.forEach(el => el.remove());
                printWindow.document.write(content.innerHTML);
                printWindow.document.write('</body></html>');
                printWindow.document.close();
                printWindow.print();
            }}
            function quickSearch() {{
                const searchTerm = document.getElementById('quick-search').value.toLowerCase();
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    const text = section.textContent.toLowerCase();
                    if (text.includes(searchTerm)) {{
                        section.style.border = '2px solid #f59e0b';
                        section.style.backgroundColor = 'rgba(245,158,11,0.1)';
                        section.scrollIntoView({{ behavior: 'smooth', block: 'nearest' }});
                    }} else {{
                        section.style.border = '';
                        section.style.backgroundColor = '';
                    }}
                }});
            }}
            function clearSearch() {{
                document.getElementById('quick-search').value = '';
                const sections = document.querySelectorAll('.card');
                sections.forEach(section => {{
                    section.style.border = '';
                    section.style.backgroundColor = '';
                }});
            }}
            let quotes = {json.dumps(self.science_quotes)};
            let currentQuote = 0;
            function rotateQuote() {{
                document.getElementById('science-quote').innerHTML = quotes[currentQuote];
                currentQuote = (currentQuote + 1) % quotes.length;
            }}
            document.addEventListener('DOMContentLoaded', function() {{
                rotateQuote();
                setInterval(rotateQuote, 10000);
            }});
        </script>
        """
        export_data_js = f"""
        <script>
            window.reportData = {{
                metadata: {{
                    genome: '{genome_name}',
                    date: '{self.metadata['analysis_date']}',
                    tool: '{self.metadata['tool_name']}',
                    version: '{self.metadata['version']}'
                }},
                summary: {json.dumps(analysis)},
                hits: {json.dumps(hits)}
            }};
        </script>
        """
        
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>EnteroMark AMRfinderPlus Analysis Report</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #92400e 0%, #d97706 50%, #f59e0b 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0, 0, 0, 0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.4);
            border: 2px solid rgba(245, 158, 11, 0.5);
        }}
        .ascii-art {{
            font-family: 'Courier New', monospace;
            font-size: 12px;
            line-height: 1.2;
            white-space: pre;
            color: #fbbf24;
            text-shadow: 0 0 10px rgba(251, 191, 36, 0.5);
            overflow-x: auto;
        }}
        .card {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.2);
        }}
        .table-container {{
            max-height: 600px;
            overflow-y: auto;
            overflow-x: auto;
            border: 1px solid #ddd;
            border-radius: 8px;
            margin: 15px 0;
        }}
        .gene-table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            min-width: 1200px;  /* Forces horizontal scroll if needed */
        }}
        .gene-table th, .gene-table td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
            font-size: 0.85rem;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            font-weight: 600;
            position: sticky;
            top: 0;
            z-index: 10;
        }}
        .class-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
        }}
        .class-table th, .class-table td {{
            padding: 15px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        .class-table th {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
        }}
        tr:hover {{ background-color: #f8f9fa; }}
        .summary-stats {{
            display: flex;
            justify-content: space-around;
            margin: 20px 0;
            flex-wrap: wrap;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #d97706 0%, #92400e 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc2626 0%, #991b1b 100%);
            color: white;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.2);
            margin: 10px;
            flex: 1;
            min-width: 200px;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin: 20px 0;
            text-align: center;
            font-style: italic;
            border-left: 4px solid #f59e0b;
        }}
        .footer {{
            background: rgba(0, 0, 0, 0.8);
            color: white;
            padding: 30px;
            border-radius: 12px;
            margin-top: 40px;
        }}
        .resistance-badge {{
            display: inline-block;
            background: #dc2626;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .critical-risk-badge {{
            display: inline-block;
            background: #8b0000;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
            font-weight: bold;
        }}
        .warning-badge {{
            display: inline-block;
            background: #f59e0b;
            color: black;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
        }}
        .present {{ background-color: #d4edda; }}
        .critical-row {{ background-color: #f8d7da; font-weight: bold; border-left: 4px solid #dc2626; }}
        .high-risk-row {{ background-color: #fff3cd; border-left: 4px solid #f59e0b; }}
        .interactive-controls {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            align-items: center;
        }}
        .search-box input {{
            width: 100%;
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
        }}
        .export-buttons button {{
            padding: 8px 16px;
            background: #d97706;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
        }}
        .export-buttons button.print {{ background: #10b981; }}
        .result-counter {{ font-size: 0.9em; color: #666; font-style: italic; }}
        .sequence-name {{ max-width: 300px; overflow-wrap: break-word; }}
        .quick-search-bar {{
            background: rgba(255, 255, 255, 0.95);
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
            display: flex;
            gap: 10px;
        }}
        .quick-search-bar input {{ flex: 1; padding: 8px 12px; border: 2px solid #d97706; border-radius: 4px; }}
        .quick-search-bar button {{ padding: 8px 20px; background: #d97706; color: white; border: none; border-radius: 4px; cursor: pointer; }}
        @media print {{
            .no-print, .interactive-controls, .quick-search-bar {{ display: none; }}
            body {{ background: white; color: black; }}
            .card {{ background: white; box-shadow: none; }}
        }}
    </style>
    {interactive_js}
    {export_data_js}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">{self.ascii_art}</div>
            </div>
            <div class="card">
                <h1 style="color: #333; font-size: 2.5em;">🧬 EnteroMark AMRfinderPlus Analysis Report</h1>
                <p style="color: #666; font-size: 1.2em;">Comprehensive Enterococcus faecium Antimicrobial Resistance Analysis</p>
                <div class="quick-search-bar no-print">
                    <input type="text" id="quick-search" placeholder="🔍 Quick search across entire report...">
                    <button onclick="quickSearch()">Search</button>
                    <button onclick="clearSearch()" style="background: #6c757d;">Clear</button>
                </div>
                <div class="export-buttons no-print" style="margin-top: 15px;">
                    <button onclick="exportToJSON('reportData', '{genome_name}_amr_report.json')">📥 Export JSON</button>
                    <button onclick="printReport()" class="print">🖨️ Print Report</button>
                </div>
            </div>
        </div>
        
        <div class="quote-container">
            <div id="science-quote" style="font-size: 1.1em;"></div>
        </div>
"""
        # Critical risk alert
        if analysis['critical_risk_genes'] > 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #dc2626; background: #fef2f2;">
            <h2 style="color: #dc2626;">🚨 CRITICAL RISK AMR GENES DETECTED</h2>
            <p><strong>{analysis['critical_risk_genes']} CRITICAL RISK antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
                <p style="color: #7f1d1d; font-weight: bold;">
                    ⚠️ These genes confer resistance to last-resort antibiotics (vancomycin, linezolid, high-level aminoglycosides) and represent a serious public health concern.
                </p>
"""
            for gene in analysis['critical_risk_list']:
                html_content += f'<span class="critical-risk-badge">🚨 {gene}</span>'
            html_content += """
            </div>
        </div>
"""
        
        html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #d97706; padding-bottom: 10px;">📊 E. faecium AMR Summary</h2>
            <div class="summary-stats">
                <div class="stat-card"><h3>Total AMR Genes</h3><p style="font-size: 2em;">{analysis['total_genes']}</p></div>
                <div class="stat-card"><h3>High Risk Genes</h3><p style="font-size: 2em;">{analysis['high_risk_genes']}</p></div>
                <div class="critical-stat-card"><h3>Critical Risk</h3><p style="font-size: 2em;">{analysis['critical_risk_genes']}</p></div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
            <p><strong>AMRfinderPlus Version:</strong> {self.metadata['amrfinder_version']}</p>
            <p><strong>Database Version:</strong> {self.metadata['database_version']}</p>
        </div>
"""
        # High-risk genes warning (non-critical)
        if analysis['high_risk_genes'] > 0 and analysis['critical_risk_genes'] == 0:
            html_content += f"""
        <div class="card" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #b45309;">⚠️ High-Risk AMR Genes Detected</h2>
            <p><strong>{analysis['high_risk_genes']} high-risk antimicrobial resistance genes found:</strong></p>
            <div style="margin: 10px 0;">
"""
            for gene in analysis['high_risk_list']:
                html_content += f'<span class="resistance-badge">{gene}</span>'
            html_content += """
            </div>
        </div>
"""
        # Resistance Mechanism Breakdown
        if any(analysis['resistance_mechanisms'].values()):
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #d97706; padding-bottom: 10px;">🔬 Resistance Mechanism Breakdown</h2>
"""
            mech = analysis['resistance_mechanisms']
            if mech['vancomycin']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #fef2f2; border-radius: 5px;">
                <strong>Vancomycin Resistance (CRITICAL):</strong> {', '.join(mech['vancomycin'])}
            </div>
"""
            if mech['linezolid']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #fef2f2; border-radius: 5px;">
                <strong>Linezolid Resistance (CRITICAL):</strong> {', '.join(mech['linezolid'])}
            </div>
"""
            if mech['aminoglycoside']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #fef2f2; border-radius: 5px;">
                <strong>High-level Aminoglycoside Resistance (CRITICAL):</strong> {', '.join(mech['aminoglycoside'])}
            </div>
"""
            if mech['macrolide_lincosamide']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #ffedd5; border-radius: 5px;">
                <strong>Macrolide/Lincosamide/Streptogramin B:</strong> {', '.join(mech['macrolide_lincosamide'])}
            </div>
"""
            if mech['tetracycline']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #e0f2fe; border-radius: 5px;">
                <strong>Tetracycline Resistance:</strong> {', '.join(mech['tetracycline'])}
            </div>
"""
            if mech['beta_lactam']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f1f5f9; border-radius: 5px;">
                <strong>Beta-lactam (including pbp5 alterations):</strong> {', '.join(mech['beta_lactam'])}
            </div>
"""
            if mech['efflux_pumps']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f3e8ff; border-radius: 5px;">
                <strong>Efflux Pumps:</strong> {', '.join(mech['efflux_pumps'])}
            </div>
"""
            if mech['other_amr']:
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 5px;">
                <strong>Other AMR Genes:</strong> {', '.join(mech['other_amr'])}
            </div>
"""
            html_content += "</div>"
        # Resistance classes summary
        if analysis['resistance_classes']:
            html_content += """
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #d97706; padding-bottom: 10px;">🧪 Resistance Classes Detected</h2>
            <div class="table-container">
                <table class="class-table">
                    <thead><tr><th>Resistance Class</th><th>Gene Count</th><th>Genes</th></tr></thead>
                    <tbody>
"""
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join(genes)
                html_content += f"""
                        <tr><td><strong>{class_name}</strong></td><td>{len(genes)}</td><td>{gene_list}</td>
                        </tr>
"""
            html_content += "</tbody></table></div></div>"
        # Detailed AMR genes table with ALL columns
        if hits:
            # Define all columns we want to display (matching the raw text file)
            all_columns = [
                'gene_symbol', 'sequence_name', 'contig_id', 'start', 'stop', 'strand',
                'class', 'subclass', 'element_type', 'element_subtype', 'scope', 'method',
                'coverage', 'identity', 'alignment_length', 'accession', 'closest_name',
                'protein_id', 'target_length', 'ref_length', 'hmm_id', 'hmm_description'
            ]
            # Generate table headers
            col_headers = {
                'gene_symbol': 'Gene Symbol',
                'sequence_name': 'Sequence Name',
                'contig_id': 'Contig',
                'start': 'Start',
                'stop': 'Stop',
                'strand': 'Strand',
                'class': 'Class',
                'subclass': 'Subclass',
                'element_type': 'Type',
                'element_subtype': 'Subtype',
                'scope': 'Scope',
                'method': 'Method',
                'coverage': 'Cov%',
                'identity': 'Id%',
                'alignment_length': 'Align Len',
                'accession': 'Accession',
                'closest_name': 'Closest Name',
                'protein_id': 'Protein ID',
                'target_length': 'Target Len',
                'ref_length': 'Ref Len',
                'hmm_id': 'HMM ID',
                'hmm_description': 'HMM Desc'
            }
            html_content += f"""
        <div class="card">
            <h2 style="color: #333; border-bottom: 2px solid #d97706; padding-bottom: 10px;">🔬 Complete AMR Data (All Fields)</h2>
            <div class="interactive-controls no-print">
                <div class="search-box"><input type="text" id="search-detailed-amr" placeholder="Search any column..." onkeyup="searchTable('detailed-amr-table', this.value)"></div>
                <div class="export-buttons"><button onclick="exportToCSV('detailed-amr-table', '{genome_name}_amr_full_data.csv')">📥 Export CSV</button></div>
                <div class="result-counter" id="result-counter-detailed-amr-table">{len(hits)} results found</div>
            </div>
            <div class="table-container">
                <table class="gene-table" id="detailed-amr-table">
                    <thead>
                        <tr>
"""
            for col in all_columns:
                html_content += f"<th>{col_headers.get(col, col)}</th>"
            html_content += """
                        </tr>
                    </thead>
                    <tbody>
"""
            for hit in hits:
                gene = hit.get('gene_symbol', '')
                row_class = "present"
                if gene in analysis['critical_risk_list']:
                    row_class = "critical-row"
                elif gene in analysis['high_risk_list']:
                    row_class = "high-risk-row"
                html_content += f'<tr class="{row_class}">'
                for col in all_columns:
                    value = hit.get(col, '')
                    # Truncate long values for display but keep full in export
                    if col in ['closest_name', 'sequence_name', 'hmm_description'] and len(value) > 10000:
                        value = value[:100000] + '...'
                    html_content += f"<td>{value}</td>"
                html_content += "</table>\n"
            html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""
        else:
            html_content += """
        <div class="card"><h2>✅ No AMR Genes Detected</h2><p>No antimicrobial resistance genes found in this Enterococcus faecium genome.</p></div>
"""
        html_content += f"""
        <div class="footer">
            <h3 style="color: #fff;">👥 Contact Information</h3>
            <p><strong>Author:</strong> Brown Beckley</p>
            <p><strong>Email:</strong> brownbeckley94@gmail.com</p>
            <p><strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub" target="_blank">https://github.com/bbeckley-hub</a></p>
            <p><strong>Affiliation:</strong> University of Ghana Medical School - Department of Medical Biochemistry</p>
            <p style="margin-top: 20px; font-size: 0.9em;"></p>
        </div>
    </div>
</body>
</html>"""
        html_file = os.path.join(output_dir, f"{genome_name}_amrfinder_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info(f"E. faecium AMRfinderPlus HTML report generated: {html_file}")
    
    def _analyze_efaecalis_amr_results(self, hits: List[Dict]) -> Dict[str, Any]:
        analysis = {
            'total_genes': len(hits),
            'resistance_classes': {},
            'total_classes': 0,
            'high_risk_genes': 0,
            'critical_risk_genes': 0,
            'high_risk_list': [],
            'critical_risk_list': [],
            'resistance_mechanisms': {
                'vancomycin': [],
                'linezolid': [],
                'aminoglycoside': [],
                'macrolide_lincosamide': [],
                'tetracycline': [],
                'beta_lactam': [],
                'efflux_pumps': [],
                'other_amr': []
            }
        }
        for hit in hits:
            gene = hit.get('gene_symbol', '')
            if gene is None or gene == '':
                continue
            gene = str(gene)
            # Categorize mechanism for E. faecium
            self._categorize_efaecalis_resistance_mechanism(gene, analysis)
            # Risk levels
            if gene in self.critical_risk_genes:
                analysis['critical_risk_genes'] += 1
                if gene not in analysis['critical_risk_list']:
                    analysis['critical_risk_list'].append(gene)
            if gene in self.high_risk_genes:
                analysis['high_risk_genes'] += 1
                if gene not in analysis['high_risk_list']:
                    analysis['high_risk_list'].append(gene)
            # Resistance class
            cls = hit.get('class', '')
            if cls:
                if cls not in analysis['resistance_classes']:
                    analysis['resistance_classes'][cls] = []
                if gene not in analysis['resistance_classes'][cls]:
                    analysis['resistance_classes'][cls].append(gene)
        analysis['total_classes'] = len(analysis['resistance_classes'])
        return analysis
    
    def _categorize_efaecalis_resistance_mechanism(self, gene: str, analysis: Dict[str, Any]):
        gene = str(gene) if gene is not None else ''
        # Vancomycin
        if any(kw in gene for kw in ('vanA', 'vanB', 'vanD', 'vanM', 'vanC', 'vanE', 'vanG', 'vanL', 'vanN')):
            analysis['resistance_mechanisms']['vancomycin'].append(gene)
        # Linezolid
        elif any(kw in gene for kw in ('optrA', 'cfr', 'poxtA')):
            analysis['resistance_mechanisms']['linezolid'].append(gene)
        # Aminoglycoside
        elif any(kw in gene for kw in ('aac', 'ant', 'aph')):
            analysis['resistance_mechanisms']['aminoglycoside'].append(gene)
        # Macrolide-lincosamide-streptogramin B
        elif any(kw in gene for kw in ('erm', 'msr', 'mef', 'mph')):
            analysis['resistance_mechanisms']['macrolide_lincosamide'].append(gene)
        # Tetracycline
        elif any(kw in gene for kw in ('tet')):
            analysis['resistance_mechanisms']['tetracycline'].append(gene)
        # Beta-lactam (pbp5)
        elif 'pbp5' in gene.lower():
            analysis['resistance_mechanisms']['beta_lactam'].append(gene)
        # Efflux pumps
        elif any(kw in gene for kw in ('msrC', 'emeA', 'efmA')):
            analysis['resistance_mechanisms']['efflux_pumps'].append(gene)
        else:
            analysis['resistance_mechanisms']['other_amr'].append(gene)
    
    def create_amr_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating E. faecium AMR summary files...")
        summary_file = os.path.join(output_base, "enteromark_amrfinder_summary.tsv")
        with open(summary_file, 'w') as f:
            f.write("Genome\tGene_Symbol\tSequence_Name\tClass\tSubclass\tCoverage\tIdentity\tScope\tElement_Type\tAccession\tContig\tStart\tStop\n")
            for genome_name, result in all_results.items():
                for hit in result['hits']:
                    row = [
                        genome_name,
                        hit.get('gene_symbol', ''),
                        hit.get('sequence_name', ''),
                        hit.get('class', ''),
                        hit.get('subclass', ''),
                        hit.get('coverage', ''),
                        hit.get('identity', ''),
                        hit.get('scope', ''),
                        hit.get('element_type', ''),
                        hit.get('accession', ''),
                        hit.get('contig_id', ''),
                        hit.get('start', ''),
                        hit.get('stop', '')
                    ]
                    f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(f"✓ Created summary: {summary_file}")
        
        stats_file = os.path.join(output_base, "enteromark_amrfinder_statistics_summary.tsv")
        with open(stats_file, 'w') as f:
            f.write("Genome\tTotal_AMR_Genes\tHigh_Risk_Genes\tCritical_Risk_Genes\tResistance_Classes\tGene_List\n")
            for genome_name, result in all_results.items():
                genes = list(set(hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol')))
                gene_list = ",".join(genes)
                high_risk_count = sum(1 for g in genes if g in self.high_risk_genes)
                critical_risk_count = sum(1 for g in genes if g in self.critical_risk_genes)
                classes = list(set(hit.get('class', '') for hit in result['hits'] if hit.get('class')))
                class_list = ",".join(classes)
                f.write(f"{genome_name}\t{result['hit_count']}\t{high_risk_count}\t{critical_risk_count}\t{class_list}\t{gene_list}\n")
        self.logger.info(f"✓ Created statistics: {stats_file}")
        
        # Create JSON master summary
        master_summary = {
            'metadata': {
                'tool': self.metadata['tool_name'],
                'version': self.metadata['version'],
                'amrfinder_version': self.metadata['amrfinder_version'],
                'database_version': self.metadata['database_version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results)
            },
            'genome_summaries': {},
            'cross_genome_patterns': {}
        }
        all_genes_by_gene = defaultdict(lambda: {'count': 0, 'genomes': set()})
        for genome_name, result in all_results.items():
            genes = [hit.get('gene_symbol', '') for hit in result['hits'] if hit.get('gene_symbol')]
            unique_genes = set(genes)
            critical_genes = [g for g in unique_genes if g in self.critical_risk_genes]
            high_risk_genes = [g for g in unique_genes if g in self.high_risk_genes and g not in self.critical_risk_genes]
            master_summary['genome_summaries'][genome_name] = {
                'total_hits': result['hit_count'],
                'unique_genes': len(unique_genes),
                'critical_genes': critical_genes,
                'high_risk_genes': high_risk_genes,
                'genes': list(unique_genes),
                'status': result['status']
            }
            for g in unique_genes:
                all_genes_by_gene[g]['count'] += 1
                all_genes_by_gene[g]['genomes'].add(genome_name)
        cross = {}
        for g, data in all_genes_by_gene.items():
            cross[g] = {
                'frequency': data['count'],
                'genomes': list(data['genomes']),
                'risk_level': 'CRITICAL' if g in self.critical_risk_genes else 'HIGH' if g in self.high_risk_genes else 'STANDARD'
            }
        master_summary['cross_genome_patterns'] = {
            'total_unique_genes': len(all_genes_by_gene),
            'genomes_with_critical': sum(1 for gs in master_summary['genome_summaries'].values() if gs['critical_genes']),
            'gene_frequency': cross
        }
        master_json = os.path.join(output_base, "enteromark_amrfinder_master_summary.json")
        with open(master_json, 'w') as f:
            json.dump(master_summary, f, indent=2)
        self.logger.info(f"✓ Created master JSON: {master_json}")
        
        # Create batch summary HTML with unique genes metrics at top and no Risk Level column
        self._create_summary_html_report(all_results, output_base)
    
    def _create_summary_html_report(self, all_results: Dict[str, Any], output_base: str):
        total_genomes = len(all_results)
        total_hits = sum(r['hit_count'] for r in all_results.values())
        
        # Compute genes per genome and gene frequency
        genes_per_genome = {}
        gene_frequency = defaultdict(set)
        for gn, res in all_results.items():
            genes = [h.get('gene_symbol', '') for h in res['hits'] if h.get('gene_symbol')]
            genes_per_genome[gn] = genes
            for g in genes:
                gene_frequency[g].add(gn)
        
        # Unique genes metrics
        all_unique_genes = set(gene_frequency.keys())
        unique_critical = {g for g in all_unique_genes if g in self.critical_risk_genes}
        unique_high = {g for g in all_unique_genes if g in self.high_risk_genes and g not in self.critical_risk_genes}
        
        # Identify critical/high‑risk genes across all genomes
        genomes_with_critical = 0
        for gn, res in all_results.items():
            genes = set(genes_per_genome[gn])
            if any(g in self.critical_risk_genes for g in genes):
                genomes_with_critical += 1
        
        # Build HTML – clean tables, no interactive controls, unique metrics at top
        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>EnteroMark AMRfinderPlus - Batch Summary</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #92400e 0%, #d97706 50%, #f59e0b 100%);
            font-family: 'Segoe UI', sans-serif;
            color: #fff;
            padding: 20px;
        }}
        .container {{ max-width: 1400px; margin: auto; }}
        .header {{ text-align: center; margin-bottom: 30px; }}
        .ascii-container {{
            background: rgba(0,0,0,0.7);
            padding: 20px;
            border-radius: 15px;
            margin-bottom: 20px;
            border: 2px solid rgba(245,158,11,0.5);
        }}
        .ascii-art {{
            font-family: monospace;
            font-size: 12px;
            white-space: pre;
            color: #fbbf24;
            text-shadow: 0 0 10px rgba(251,191,36,0.5);
            overflow-x: auto;
        }}
        .card {{
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 25px;
            margin: 20px 0;
            border-radius: 12px;
        }}
        .card h2 {{
            color: #b45309;
            border-bottom: 3px solid #d97706;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        .summary-stats {{
            display: flex;
            gap: 20px;
            justify-content: center;
            flex-wrap: wrap;
            margin: 20px 0;
        }}
        .stat-card {{
            background: linear-gradient(135deg, #d97706 0%, #92400e 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            flex: 1;
            min-width: 200px;
        }}
        .critical-stat-card {{
            background: linear-gradient(135deg, #dc2626 0%, #991b1b 100%);
        }}
        .table-container {{
            max-height: 600px;
            overflow-y: auto;
            overflow-x: auto;
            margin-top: 20px;
            border: 1px solid #ddd;
            border-radius: 8px;
        }}
        .gene-table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            min-width: 800px;
        }}
        .gene-table th, .gene-table td {{
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e5e7eb;
        }}
        .gene-table th {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            position: sticky;
            top: 0;
        }}
        .gene-table tr:nth-child(even) {{ background-color: #f8fafc; }}
        .gene-table tr:hover {{ background-color: #e0f2fe; }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .high-risk {{ background-color: #fef3c7; }}
        .risk-badge {{
            display: inline-block;
            background: #dc2626;
            color: white;
            padding: 4px 8px;
            border-radius: 12px;
            margin: 2px;
            font-size: 0.8em;
        }}
        .warning-badge {{
            background: #f59e0b;
            color: black;
        }}
        .footer {{
            background: rgba(0,0,0,0.8);
            color: white;
            padding: 20px;
            border-radius: 12px;
            margin-top: 40px;
            text-align: center;
        }}
        .timestamp {{ color: #fbbf24; }}
        .authorship {{ margin-top: 15px; font-size: 0.9em; }}
    </style>
</head>
<body>
<div class="container">
    <div class="header">
        <div class="ascii-container">
            <div class="ascii-art">{self.ascii_art}</div>
        </div>
        <div class="card">
            <h1>🧬 EnteroMark AMRfinderPlus – Batch Summary</h1>
            <p>Enterococcus faecium Antimicrobial Resistance Analysis</p>
            <p><strong>Date:</strong> {self.metadata['analysis_date']} | <strong>Tool:</strong> {self.metadata['tool_name']} v{self.metadata['version']}</p>
        </div>
    </div>
    
    <div class="card">
        <h2>📊 Overall Summary</h2>
        <div class="summary-stats">
            <div class="stat-card"><div class="stat-label">Total Genomes</div><div class="stat-value" style="font-size:2em;">{total_genomes}</div></div>
            <div class="stat-card"><div class="stat-label">Total AMR Genes</div><div class="stat-value" style="font-size:2em;">{total_hits}</div></div>
            <div class="stat-card critical-stat-card"><div class="stat-label">Genomes with Critical Risk</div><div class="stat-value" style="font-size:2em;">{genomes_with_critical}</div></div>
        </div>
    </div>
    
    <div class="card">
        <h2>📊 Unique Genes Across All Genomes</h2>
        <div class="summary-stats">
            <div class="stat-card">Total unique genes: {len(all_unique_genes)}</div>
            <div class="stat-card">Unique critical genes: {len(unique_critical)}</div>
            <div class="stat-card">Unique high-risk genes: {len(unique_high)}</div>
        </div>
        <div>
            <p><strong>Critical genes detected:</strong> {', '.join(sorted(unique_critical)) if unique_critical else 'None'}</p>
            <p><strong>High-risk genes detected (non-critical):</strong> {', '.join(sorted(unique_high)) if unique_high else 'None'}</p>
        </div>
    </div>
"""
        if unique_critical:
            html += f"""
    <div class="card" style="border-left: 4px solid #dc2626;">
        <h2 style="color: #dc2626;">🚨 CRITICAL RISK GENES ACROSS ALL GENOMES</h2>
        <div>
            <p><strong>{len(unique_critical)} unique critical risk genes found in {genomes_with_critical} genomes:</strong></p>
            <div>"""
            for g in sorted(unique_critical):
                html += f'<span class="risk-badge">🚨 {g}</span>'
            html += "</div></div></div>"
        
        # Genes by Genome table
        html += """
    <div class="card">
        <h2>🔍 Genes by Genome</h2>
        <div class="table-container">
            <table class="gene-table">
                <thead><tr><th>Genome</th><th>Gene Count</th><th>Genes Detected</th></tr></thead>
                <tbody>
"""
        for gn in sorted(genes_per_genome.keys()):
            genes = genes_per_genome[gn]
            row_class = "critical" if any(g in self.critical_risk_genes for g in genes) else "high-risk" if any(g in self.high_risk_genes for g in genes) else ""
            html += f"<tr class='{row_class}'> <td><strong>{gn}</strong></td><td>{len(genes)}</td><td>{', '.join(sorted(genes))}</td></tr>"
        html += """
                </tbody>
            </table>
        </div>
    </div>
    
    <div class="card">
        <h2>📈 Gene Frequency</h2>
        <div class="table-container">
            <table class="gene-table">
                <thead> <tr><th>Gene</th><th>Frequency</th><th>Prevalence</th><th>Genomes</th> </tr> </thead>
                <tbody>
"""
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            freq = len(genomes)
            pct = (freq / total_genomes) * 100
            if pct >= 75:
                prev = '<span class="risk-badge">Very High</span>'
            elif pct >= 50:
                prev = '<span class="risk-badge warning-badge">High</span>'
            elif pct >= 25:
                prev = '<span class="risk-badge warning-badge">Medium</span>'
            elif pct >= 10:
                prev = '<span class="risk-badge">Low</span>'
            else:
                prev = '<span class="risk-badge">Rare</span>'
            html += f"<tr> <td><strong>{gene}</strong></td><td>{freq} ({pct:.1f}%)</td><td>{prev}</td><td>{', '.join(sorted(genomes))}</td></tr>"
        html += f"""
                </tbody>
            </table>
        </div>
    </div>
    
    <div class="footer">
        <p><strong>ENTEROMARK</strong> – AMRfinderPlus Batch Analysis for Enterococcus faecium</p>
        <p class="timestamp">Generated: {self.metadata['analysis_date']}</p>
        <div class="authorship">
            <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
            <p>Email: brownbeckley94@gmail.com</p>
            <p>University of Ghana Medical School – Department of Medical Biochemistry</p>
        </div>
    </div>
</div>
</body>
</html>"""
        html_file = os.path.join(output_base, "enteromark_amrfinder_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html)
        self.logger.info(f"✓ Created summary HTML report: {html_file}")
    
    def process_single_genome(self, genome_file: str, output_base: str = "enteromark_amr_results") -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        os.makedirs(results_dir, exist_ok=True)
        self.logger.info(f"=== PROCESSING GENOME: {genome_name} ===")
        result = self.run_amrfinder_single_genome(genome_file, results_dir)
        if result['status'] == 'success':
            try:
                self._create_amrfinder_html_report(genome_name, result['hits'], results_dir)
            except Exception as e:
                self.logger.error(f"Failed to generate HTML report for {genome_name}: {e}")
        status = "✓" if result['status'] == 'success' else "✗"
        self.logger.info(f"{status} {genome_name}: {result['hit_count']} AMR hits")
        return result
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "enteromark_amr_results") -> Dict[str, Any]:
        if not self.check_amrfinder_installed():
            raise RuntimeError("Bundled AMRfinderPlus not properly installed")
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa",
                          f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        genome_files = []
        for pat in fasta_patterns:
            genome_files.extend(glob.glob(pat))
        genome_files = list(set(genome_files))
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching: {genome_pattern}")
        self.logger.info(f"Found {len(genome_files)} genomes")
        os.makedirs(output_base, exist_ok=True)
        
        # RAM-aware concurrency
        max_concurrent = self._calculate_max_concurrent(len(genome_files))
        
        self.logger.info(f"🚀 MAXIMUM SPEED: Using {max_concurrent} concurrent jobs")
        with ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            future_to_genome = {executor.submit(self.process_single_genome, g, output_base): g for g in genome_files}
            all_results = {}
            for future in as_completed(future_to_genome):
                genome = future_to_genome[future]
                try:
                    result = future.result()
                    all_results[result['genome']] = result
                    self.logger.info(f"✓ COMPLETED: {result['genome']} ({result['hit_count']} hits)")
                except Exception as e:
                    self.logger.error(f"✗ FAILED: {genome} - {e}")
                    all_results[Path(genome).stem] = {'genome': Path(genome).stem, 'hits': [], 'hit_count': 0, 'status': 'failed'}
        self.create_amr_summary(all_results, output_base)
        self.logger.info("=== E. FAECIUM AMR ANALYSIS COMPLETE ===")
        return all_results

def main():
    parser = argparse.ArgumentParser(
        description='EnteroMark AMRfinderPlus - Enterococcus faecium Antimicrobial Resistance - MAXIMUM SPEED',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python enteromark_amr.py "*.fasta"
  python enteromark_amr.py "*.fna" --output results --cpus 8
  python enteromark_amr.py --update-db
  python enteromark_amr.py --db-version
        """
    )
    parser.add_argument('pattern', nargs='?', help='File pattern for genomes (e.g., "*.fasta")')
    parser.add_argument('--cpus', '-c', type=int, default=None, help='CPU cores (auto-detect if not set)')
    parser.add_argument('--output', '-o', default='enteromark_amr_results', help='Output directory')
    parser.add_argument('--update-db', action='store_true', help='Update AMRfinderPlus database to latest version and exit')
    parser.add_argument('--db-version', action='store_true', help='Show current database version and exit')
    
    args = parser.parse_args()
    
    # Handle database operations without requiring pattern
    executor = EnteroMarkAMRfinder(cpus=args.cpus)
    if args.update_db:
        print("Updating AMRfinderPlus database...")
        success = executor.update_database()
        if success:
            print("Database updated successfully.")
        else:
            print("Database update failed.")
        sys.exit(0)
    if args.db_version:
        print(executor.metadata['database_version'])   
        sys.exit(0)
    
    # For analysis, pattern is required
    if not args.pattern:
        parser.error("Please provide a file pattern for genomes (or use --update-db / --db-version)")
    
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        executor.logger.info("\n" + "="*50)
        executor.logger.info("🧫 EnteroMark AMRfinderPlus FINAL SUMMARY")
        executor.logger.info("="*50)
        total_hits = sum(r['hit_count'] for r in results.values())
        total_crit = 0
        total_high = 0
        for r in results.values():
            genes = [h.get('gene_symbol') for h in r['hits'] if h.get('gene_symbol')]
            total_crit += sum(1 for g in genes if g in executor.critical_risk_genes)
            total_high += sum(1 for g in genes if g in executor.high_risk_genes)
        executor.logger.info(f"Genomes processed: {len(results)}")
        executor.logger.info(f"Total AMR hits: {total_hits}")
        executor.logger.info(f"High‑risk genes detected: {total_high}")
        executor.logger.info(f"CRITICAL RISK genes detected: {total_crit}")
        executor.logger.info(f"Results saved to: {args.output}")
        executor.logger.info(f"CPU cores used: {executor.cpus}")
        executor.logger.info(f"Bundled AMRfinderPlus: {executor.metadata['amrfinder_version']}")
        executor.logger.info(f"Bundled database: {executor.metadata['database_version']}")
    except Exception as e:
        executor.logger.error(f"Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()