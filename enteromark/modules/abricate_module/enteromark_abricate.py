#!/usr/bin/env python3
"""
EnteroMark ABRicate Standalone Module
Comprehensive ABRicate analysis for Enterococcus faecium with HTML, TSV, and JSON reporting
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Date: 2026-04-12
Version: 1.0.0 (E. faecium-optimized, complete)
"""

import subprocess
import sys
import os
import glob
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Dict, Any, Set, Tuple
import argparse
import re
from datetime import datetime
import psutil
import json
import random
from collections import defaultdict, Counter

class EnteroMarkAbricateExecutor:
    """ABRicate executor for Enterococcus faecium with comprehensive reporting - MAXIMUM SPEED"""
    
    def __init__(self, cpus: int = None):
        self.logger = self._setup_logging()
        self.available_ram = self._get_available_ram()
        self.cpus = self._calculate_optimal_cpus(cpus)
        
        # Databases to be used
        self.required_databases = [
            'ncbi', 'card', 'resfinder', 'vfdb', 'argannot',
            'plasmidfinder', 'megares', 'bacmet2'
        ]
        
        # ==================== COMPREHENSIVE E. FAECIUM GENE DICTIONARIES ====================
        self.critical_resistance_genes = {
            'vanA', 'vanB', 'vanD', 'vanM',
            'optrA', 'cfr', 'cfrB', 'poxtA',
            'aac(6\')-Ii', 'ant(6)-Ia', 'aph(3\')-III', 'aac(6\')-Ie-aph(2\'\')-Ia'
        }
        
        self.high_risk_virulence_genes = {
            'esp', 'ace', 'asal', 'gelE', 'cylA', 'cylB', 'cylL', 'cylM', 'cylR',
            'efaA', 'fsrA', 'fsrB', 'fsrC', 'hyl', 'scm', 'acm', 'ecbA', 'ecbB'
        }
        
        self.beta_lactamase_genes = {
            'pbp5', 'pbp5_mutation', 'pbp5_S462A', 'pbp5_E629V', 'pbp5_M485A', 'pbp5_M426I'
        }
        
        self.metadata = {
            "tool_name": "EnteroMark ABRicate",
            "version": "1.0.0",
            "authors": ["Brown Beckley"],
            "email": "brownbeckley94@gmail.com",
            "github": "https://github.com/bbeckley-hub",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        }
        
        self.science_quotes = [
            {"text": "The important thing is not to stop questioning. Curiosity has its own reason for existing.", "author": "Albert Einstein"},
            {"text": "Science is not only a disciple of reason but also one of romance and passion.", "author": "Stephen Hawking"},
            {"text": "Somewhere, something incredible is waiting to be known.", "author": "Carl Sagan"},
            {"text": "In science, there are no shortcuts to truth.", "author": "Karl Popper"},
            {"text": "Science knows no country, because knowledge belongs to humanity.", "author": "Louis Pasteur"},
            {"text": "The science of today is the technology of tomorrow.", "author": "Edward Teller"},
            {"text": "Research is what I'm doing when I don't know what I'm doing.", "author": "Wernher von Braun"},
            {"text": "EnteroMark turns genomic complexity into actionable insights for VRE surveillance.", "author": "Brown Beckley"}
        ]
    
    def get_random_quote(self):
        return random.choice(self.science_quotes)
    
    def _setup_logging(self):
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        return logging.getLogger(__name__)
    
    def _get_available_ram(self) -> int:
        try:
            ram_gb = psutil.virtual_memory().available / (1024 ** 3)
            return ram_gb
        except Exception:
            return 8
    
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
    
    def check_abricate_installed(self) -> bool:
        try:
            result = subprocess.run(['abricate', '--version'],
                                    capture_output=True, text=True, check=True)
            version_line = result.stdout.strip()
            self.logger.info("ABRicate version: %s", version_line)
            version_match = re.search(r'(\d+\.\d+\.\d+)', version_line)
            if version_match and version_match.group(1) >= "1.2.0":
                self.logger.info("✓ ABRicate version meets requirement (>=1.2.0)")
                return True
            self.logger.info("✓ ABRicate installed")
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.logger.error("ABRicate not found. Please install with: conda install -c bioconda abricate")
            return False
    
    def setup_abricate_databases(self):
        self.logger.info("Setting up ABRicate databases...")
        available_dbs = []
        missing_dbs = []
        try:
            check_result = subprocess.run(['abricate', '--list'],
                                          capture_output=True, text=True, check=True)
            for db in self.required_databases:
                if db in check_result.stdout:
                    self.logger.info("✓ Database available: %s", db)
                    available_dbs.append(db)
                else:
                    self.logger.warning("Database not available: %s", db)
                    missing_dbs.append(db)
            for db in missing_dbs:
                self.logger.info("Attempting to setup database: %s", db)
                try:
                    subprocess.run(['abricate', '--setupdb', '--db', db],
                                   capture_output=True, text=True, check=True)
                    self.logger.info("✓ Database setup completed: %s", db)
                    available_dbs.append(db)
                except subprocess.CalledProcessError as e:
                    self.logger.error("Failed to setup database %s: %s", db, e.stderr)
            self.required_databases = available_dbs
            self.logger.info("Using databases: %s", ", ".join(self.required_databases))
        except Exception as e:
            self.logger.error("Error setting up databases: %s", e)
    
    def run_abricate_single_db(self, genome_file: str, database: str, output_dir: str) -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        output_file = os.path.join(output_dir, f"abricate_{database}.txt")
        cmd = [
            'abricate',
            genome_file,
            '--db', database,
            '--minid', '80',
            '--mincov', '80'
        ]
        self.logger.info("Running ABRicate: %s --db %s", genome_name, database)
        try:
            with open(output_file, 'w') as outfile:
                subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True, check=True)
            hits = self._parse_abricate_output(output_file)
            self._create_database_html_report(genome_name, database, hits, output_dir)
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': hits,
                'hit_count': len(hits),
                'status': 'success'
            }
        except subprocess.CalledProcessError as e:
            self.logger.error("ABRicate failed for %s on %s: %s", database, genome_name, e.stderr)
            return {
                'database': database,
                'genome': genome_name,
                'output_file': output_file,
                'hits': [],
                'hit_count': 0,
                'status': 'failed'
            }
    
    def _parse_abricate_output(self, abricate_file: str) -> List[Dict]:
        hits = []
        try:
            with open(abricate_file, 'r') as f:
                lines = f.readlines()
            if not lines:
                return hits
            headers = []
            data_lines = []
            for line in lines:
                if line.startswith('#FILE') and not headers:
                    headers = line.strip().replace('#', '').split('\t')
                elif line.strip() and not line.startswith('#'):
                    data_lines.append(line.strip())
            if not headers:
                return hits
            expected_columns = len(headers)
            for line_num, line in enumerate(data_lines, 1):
                parts = line.split('\t')
                if len(parts) > expected_columns:
                    combined_parts = parts[:expected_columns-1]
                    combined_parts.append('\t'.join(parts[expected_columns-1:]))
                    parts = combined_parts
                elif len(parts) < expected_columns:
                    parts.extend([''] * (expected_columns - len(parts)))
                if len(parts) == expected_columns:
                    hit = {}
                    for i, header in enumerate(headers):
                        hit[header] = parts[i] if i < len(parts) else ''
                    processed_hit = {
                        'file': hit.get('FILE', ''),
                        'sequence': hit.get('SEQUENCE', ''),
                        'start': hit.get('START', ''),
                        'end': hit.get('END', ''),
                        'strand': hit.get('STRAND', ''),
                        'gene': hit.get('GENE', ''),
                        'coverage': hit.get('COVERAGE', ''),
                        'coverage_map': hit.get('COVERAGE_MAP', ''),
                        'gaps': hit.get('GAPS', ''),
                        'coverage_percent': hit.get('%COVERAGE', ''),
                        'identity_percent': hit.get('%IDENTITY', ''),
                        'database': hit.get('DATABASE', ''),
                        'accession': hit.get('ACCESSION', ''),
                        'product': hit.get('PRODUCT', ''),
                        'resistance': hit.get('RESISTANCE', '')
                    }
                    hits.append(processed_hit)
                else:
                    self.logger.warning("Line %d has %d parts, expected %d", line_num, len(parts), expected_columns)
        except Exception as e:
            self.logger.error("Error parsing %s: %s", abricate_file, e)
        self.logger.info("Parsed %d hits from %s", len(hits), abricate_file)
        return hits
    
    def _create_database_html_report(self, genome_name: str, database: str, hits: List[Dict], output_dir: str):
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        table_rows = ''
        for hit in hits:
            gene = hit['gene']
            row_class = "present"
            if any(crit in gene for crit in self.critical_resistance_genes):
                row_class = "critical"
            elif any(vir in gene for vir in self.high_risk_virulence_genes):
                row_class = "high-risk"
            table_rows += f"""
                    <tr class="{row_class}">
                        <td>{hit['file']}</td>
                        <td>{hit['sequence']}</td>
                        <td>{hit['start']}</td>
                        <td>{hit['end']}</td>
                        <td>{hit['strand']}</td>
                        <td><strong>{hit['gene']}</strong></td>
                        <td>{hit['coverage']}</td>
                        <td class="mono">{hit['coverage_map']}</td>
                        <td>{hit['gaps']}</td>
                        <td>{hit['coverage_percent']}%</td>
                        <td>{hit['identity_percent']}%</td>
                        <td>{hit['database']}</td>
                        <td>{hit['accession']}</td>
                        <td class="product-cell">{hit['product']}</td>
                        <td class="resistance-cell">{hit['resistance']}</td>
                    </tr>
"""
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroMark - ABRicate {database.upper()} Database Report</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            background: linear-gradient(135deg, #92400e 0%, #d97706 50%, #f59e0b 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
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
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #fbbf24;
            text-shadow: 0 0 10px rgba(251, 191, 36, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
        }}
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #b45309;
            border-bottom: 3px solid #d97706;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }}
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
        }}
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 20px 0;
        }}
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 12px;
            min-width: 1400px;
        }}
        .data-table th {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
            position: sticky;
            top: 0;
        }}
        .data-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
            vertical-align: top;
        }}
        .data-table tr:nth-child(even) {{
            background-color: #fef3c7;
        }}
        .data-table tr:hover {{
            background-color: #fffbeb;
        }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .high-risk {{ background-color: #fef3c7; }}
        .present {{ background-color: #f8fafc; }}
        .mono {{
            font-family: 'Courier New', monospace;
            font-size: 11px;
        }}
        .product-cell, .resistance-cell {{
            white-space: normal;
            word-wrap: break-word;
            max-width: 300px;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .metrics-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">
╔═══════════════════════════════════════════════════════════════════════════════╗
║   ███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███╗   ███╗ █████╗ ██████╗ ██╗  ██╗
║   ██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗████╗ ████║██╔══██╗██╔══██╗██║ ██╔╝
║   █████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║██╔████╔██║███████║██████╔╝█████╔╝ 
║   ██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║██║╚██╔╝██║██╔══██║██╔══██╗██╔═██╗ 
║   ███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝██║ ╚═╝ ██║██║  ██║██║  ██║██║  ██╗
║   ╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝
║                     EnteroMark - E. faecium ABRicate Analysis                 ║
╚═══════════════════════════════════════════════════════════════════════════════╝
                </div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Database Information</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Database Name</div>
                    <div class="metric-value">{database.upper()}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Genome</div>
                    <div class="metric-value">{genome_name}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Analysis Date</div>
                    <div class="metric-value">{datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Total Hits</div>
                    <div class="metric-value">{len(hits)}</div>
                </div>
            </div>
        </div>
"""
        if hits:
            html_content += """
        <div class="report-section">
            <h2>🔍 Complete ABRicate Results (All Columns)</h2>
            <div class="table-responsive">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>FILE</th>
                            <th>SEQUENCE</th>
                            <th>START</th>
                            <th>END</th>
                            <th>STRAND</th>
                            <th>GENE</th>
                            <th>COVERAGE</th>
                            <th>COVERAGE_MAP</th>
                            <th>GAPS</th>
                            <th>%COVERAGE</th>
                            <th>%IDENTITY</th>
                            <th>DATABASE</th>
                            <th>ACCESSION</th>
                            <th>PRODUCT</th>
                            <th>RESISTANCE</th>
                        </tr>
                    </thead>
                    <tbody>
"""
            html_content += table_rows
            html_content += """
                    </tbody>
                </table>
            </div>
        </div>
"""
        else:
            html_content += f"""
        <div class="report-section">
            <h2>✅ No Genes Detected</h2>
            <p>No significant hits found in the {database.upper()} database.</p>
        </div>
"""
        html_content += f"""
        <div class="footer">
            <p><strong>ENTEROMARK</strong> - ABRicate Analysis Module</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
    <script>
        const quotes = {json.dumps(self.science_quotes)};
        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');
        function getRandomQuote() {{ return quotes[Math.floor(Math.random() * quotes.length)]; }}
        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>"""
        html_file = os.path.join(output_dir, f"abricate_{database}_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Individual database report: %s", html_file)
    
    def analyze_efaecalis_genes(self, all_hits: List[Dict]) -> Dict[str, Any]:
        """Classify hits into critical resistance, high-risk virulence, and other."""
        analysis = {
            'critical_resistance_genes': [],
            'high_risk_virulence_genes': [],
            'beta_lactamase_genes': [],
            'other_genes': [],
            'resistance_classes': {},
            'total_critical_resistance': 0,
            'total_high_risk_virulence': 0,
            'total_beta_lactamase': 0,
            'total_other': 0,
            'total_hits': len(all_hits)
        }
        for hit in all_hits:
            gene = hit['gene']
            if any(crit in gene for crit in self.critical_resistance_genes):
                analysis['critical_resistance_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'CRITICAL-RES'
                })
            elif any(vir in gene for vir in self.high_risk_virulence_genes):
                analysis['high_risk_virulence_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'HIGH-VIRULENCE'
                })
            elif any(bla in gene for bla in self.beta_lactamase_genes):
                analysis['beta_lactamase_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'risk_level': 'BETA-LACTAMASE'
                })
            else:
                analysis['other_genes'].append({
                    'gene': gene,
                    'product': hit['product'],
                    'database': hit['database'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent']
                })
            res_class = self._classify_resistance(hit['product'])
            if res_class:
                if res_class not in analysis['resistance_classes']:
                    analysis['resistance_classes'][res_class] = []
                if gene not in [g['gene'] for g in analysis['resistance_classes'][res_class]]:
                    analysis['resistance_classes'][res_class].append({'gene': gene, 'product': hit['product']})
        analysis['total_critical_resistance'] = len(analysis['critical_resistance_genes'])
        analysis['total_high_risk_virulence'] = len(analysis['high_risk_virulence_genes'])
        analysis['total_beta_lactamase'] = len(analysis['beta_lactamase_genes'])
        analysis['total_other'] = len(analysis['other_genes'])
        return analysis
    
    def _classify_resistance(self, product: str) -> str:
        product_lower = product.lower()
        if any(term in product_lower for term in ['vancomycin', 'van', 'glycopeptide']):
            return 'Vancomycin resistance'
        elif any(term in product_lower for term in ['linezolid', 'oxazolidinone', 'optr', 'cfr', 'poxt']):
            return 'Linezolid resistance'
        elif any(term in product_lower for term in ['aminoglycoside', 'aac', 'ant', 'aph']):
            return 'Aminoglycoside resistance'
        elif any(term in product_lower for term in ['beta-lactam', 'penicillin', 'pbp5']):
            return 'Beta-lactam resistance'
        elif any(term in product_lower for term in ['macrolide', 'erm', 'msr', 'mef']):
            return 'Macrolide resistance'
        elif any(term in product_lower for term in ['tetracycline', 'tet']):
            return 'Tetracycline resistance'
        elif any(term in product_lower for term in ['fluoroquinolone', 'qnr']):
            return 'Fluoroquinolone resistance'
        else:
            return 'Other resistance'
    
    def create_comprehensive_html_report(self, genome_name: str, results: Dict, output_dir: str):
        all_hits = []
        for db_result in results.values():
            all_hits.extend(db_result['hits'])
        analysis = self.analyze_efaecalis_genes(all_hits)
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
        table_rows = ''
        for hit in all_hits:
            gene = hit['gene']
            row_class = "present"
            if any(crit in gene for crit in self.critical_resistance_genes):
                row_class = "critical"
            elif any(vir in gene for vir in self.high_risk_virulence_genes):
                row_class = "high-risk"
            table_rows += f"""
                    <tr class="{row_class}">
                        <td>{hit['file']}</td>
                        <td>{hit['sequence']}</td>
                        <td>{hit['start']}</td>
                        <td>{hit['end']}</td>
                        <td>{hit['strand']}</td>
                        <td><strong>{hit['gene']}</strong></td>
                        <td>{hit['coverage']}</td>
                        <td class="mono">{hit['coverage_map']}</td>
                        <td>{hit['gaps']}</td>
                        <td>{hit['coverage_percent']}%</td>
                        <td>{hit['identity_percent']}%</td>
                        <td>{hit['database']}</td>
                        <td>{hit['accession']}</td>
                        <td class="product-cell">{hit['product']}</td>
                        <td class="resistance-cell">{hit['resistance']}</td>
                    </tr>
"""
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroMark - ABRicate Analysis Report</title>
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
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #fbbf24;
            text-shadow: 0 0 10px rgba(251, 191, 36, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{ font-size: 18px; font-style: italic; margin-bottom: 10px; }}
        .quote-author {{ font-size: 14px; color: #fbbf24; font-weight: bold; }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #b45309;
            border-bottom: 3px solid #d97706;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .metric-label {{ font-size: 14px; opacity: 0.9; margin-bottom: 5px; }}
        .metric-value {{ font-size: 24px; font-weight: bold; }}
        .risk-badge {{
            display: inline-block;
            background: #dc2626;
            color: white;
            padding: 5px 10px;
            border-radius: 15px;
            margin: 2px;
            font-size: 0.9em;
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
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 20px 0;
        }}
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 12px;
            min-width: 1400px;
        }}
        .data-table th {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
            position: sticky;
            top: 0;
        }}
        .data-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
            vertical-align: top;
        }}
        .data-table tr:nth-child(even) {{ background-color: #fef3c7; }}
        .data-table tr:hover {{ background-color: #fffbeb; }}
        .critical {{ background-color: #fee2e2; font-weight: bold; }}
        .high-risk {{ background-color: #fef3c7; }}
        .present {{ background-color: #f8fafc; }}
        .mono {{ font-family: 'Courier New', monospace; font-size: 11px; }}
        .product-cell, .resistance-cell {{ white-space: normal; word-wrap: break-word; max-width: 300px; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{ color: #fbbf24; font-weight: bold; }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .metrics-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">
╔═══════════════════════════════════════════════════════════════════════════════╗
║   ███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███╗   ███╗ █████╗ ██████╗ ██╗  ██╗
║   ██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗████╗ ████║██╔══██╗██╔══██╗██║ ██╔╝
║   █████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║██╔████╔██║███████║██████╔╝█████╔╝ 
║   ██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║██║╚██╔╝██║██╔══██║██╔══██╗██╔═██╗ 
║   ███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝██║ ╚═╝ ██║██║  ██║██║  ██║██║  ██╗
║   ╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝
║                     EnteroMark - E. faecium ABRicate Analysis                 ║
╚═══════════════════════════════════════════════════════════════════════════════╝
                </div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 E. faecium AMR/Virulence Summary</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-label">Total Genes</div>
                    <div class="metric-value">{analysis['total_hits']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Critical Resistance</div>
                    <div class="metric-value">{analysis['total_critical_resistance']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">High-Risk Virulence</div>
                    <div class="metric-value">{analysis['total_high_risk_virulence']}</div>
                </div>
                <div class="metric-card">
                    <div class="metric-label">Beta-Lactamases</div>
                    <div class="metric-value">{analysis['total_beta_lactamase']}</div>
                </div>
            </div>
            <p><strong>Genome:</strong> {genome_name}</p>
            <p><strong>Tool Version:</strong> {self.metadata['version']}</p>
        </div>
"""
        if analysis['critical_resistance_genes']:
            html_content += """
        <div class="report-section" style="border-left: 4px solid #dc2626;">
            <h2 style="color: #dc2626;">⚠️ CRITICAL RESISTANCE GENES DETECTED</h2>
            <div style="margin: 10px 0;">
"""
            for gene_info in analysis['critical_resistance_genes']:
                html_content += f'<span class="risk-badge">{gene_info["gene"]}</span> '
            html_content += """
            </div>
        </div>
"""
        if analysis['high_risk_virulence_genes']:
            html_content += """
        <div class="report-section" style="border-left: 4px solid #f59e0b;">
            <h2 style="color: #f59e0b;">🟡 HIGH-RISK VIRULENCE GENES DETECTED</h2>
            <div style="margin: 10px 0;">
"""
            for gene_info in analysis['high_risk_virulence_genes']:
                html_content += f'<span class="warning-badge">{gene_info["gene"]}</span> '
            html_content += """
            </div>
        </div>
"""
        if analysis['resistance_classes']:
            html_content += """
        <div class="report-section">
            <h2>🧪 Resistance Classes Detected</h2>
"""
            for class_name, genes in analysis['resistance_classes'].items():
                gene_list = ", ".join([g['gene'] for g in genes])
                html_content += f"""
            <div style="margin: 10px 0; padding: 10px; background: #f8f9fa; border-radius: 8px;">
                <strong style="color: #667eea;">{class_name}</strong> ({len(genes)} genes)
                <br><span style="color: #666; font-size: 0.9em;">{gene_list}</span>
            </div>
"""
            html_content += "</div>"
        
        if analysis['critical_resistance_genes']:
            html_content += """
        <div class="report-section">
            <h2>🔴 Critical Resistance Genes</h2>
            <div class="table-responsive">
                <table class="data-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['critical_resistance_genes']:
                html_content += f"""
                    <tr class="critical">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            html_content += "</tbody></table></div></div>"
        
        if analysis['high_risk_virulence_genes']:
            html_content += """
        <div class="report-section">
            <h2>🟡 High-Risk Virulence Genes</h2>
            <div class="table-responsive">
                <table class="data-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['high_risk_virulence_genes']:
                html_content += f"""
                    <tr class="high-risk">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            html_content += "</tbody></table></div></div>"
        
        if analysis['beta_lactamase_genes']:
            html_content += """
        <div class="report-section">
            <h2>🔵 Beta-Lactamase Genes (non-critical)</h2>
            <div class="table-responsive">
                <table class="data-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['beta_lactamase_genes'][:1000]:
                html_content += f"""
                    <tr class="present">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            if len(analysis['beta_lactamase_genes']) > 1000:
                html_content += f"<tr><td colspan='5'>... and {len(analysis['beta_lactamase_genes'])-1000} more</td></tr>"
            html_content += "</tbody></table></div></div>"
        
        if analysis['other_genes']:
            html_content += """
        <div class="report-section">
            <h2>🔵 Other Detected Genes</h2>
            <div class="table-responsive">
                <table class="data-table">
                    <thead><tr><th>Gene</th><th>Product</th><th>Database</th><th>Coverage</th><th>Identity</th></tr></thead>
                    <tbody>
"""
            for g in analysis['other_genes'][:5000]:
                html_content += f"""
                    <tr class="present">
                        <td><strong>{g['gene']}</strong></td>
                        <td class="product-cell">{g['product']}</td>
                        <td>{g['database']}</td>
                        <td>{g['coverage']}%</td>
                        <td>{g['identity']}%</td>
                    </tr>
"""
            if len(analysis['other_genes']) > 5000:
                html_content += f"<tr><td colspan='5'>... and {len(analysis['other_genes'])-5000} more</td></tr>"
            html_content += "</tbody><tr></div></div>"
        
        html_content += """
        <div class="report-section">
            <h2>📋 Complete ABRicate Data (All Columns)</h2>
            <div class="table-responsive">
                <table class="data-table">
                    <thead>
                        <tr>
                            <th>FILE</th>
                            <th>SEQUENCE</th>
                            <th>START</th>
                            <th>END</th>
                            <th>STRAND</th>
                            <th>GENE</th>
                            <th>COVERAGE</th>
                            <th>COVERAGE_MAP</th>
                            <th>GAPS</th>
                            <th>%COVERAGE</th>
                            <th>%IDENTITY</th>
                            <th>DATABASE</th>
                            <th>ACCESSION</th>
                            <th>PRODUCT</th>
                            <th>RESISTANCE</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        html_content += table_rows
        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🗃️ Database Results Summary</h2>
            <div class="table-responsive">
                <table class="data-table">
                    <thead><tr><th>Database</th><th>Hits</th><th>Status</th></tr></thead>
                    <tbody>
"""
        for db, result in results.items():
            status_icon = "✅" if result['status'] == 'success' else "❌"
            html_content += f"""
                    <tr>
                        <td>{db}</td>
                        <td>{result['hit_count']}</td>
                        <td>{status_icon} {result['status']}</td>
                    </tr>
"""
        html_content += f"""
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>ENTEROMARK</strong> - ABRicate Analysis Module</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p><strong>Technical Support & Inquiries:</strong></p>
                <p>Author: Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>Affiliation: University of Ghana Medical School - Department of Medical Biochemistry</p>
            </div>
        </div>
    </div>
    <script>
        const quotes = {json.dumps(self.science_quotes)};
        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');
        function getRandomQuote() {{ return quotes[Math.floor(Math.random() * quotes.length)]; }}
        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>"""
        html_file = os.path.join(output_dir, f"{genome_name}_comprehensive_abricate_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Comprehensive E. faecium HTML report generated: %s", html_file)
    
    def create_database_summaries(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating database summary files and HTML reports...")
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = []
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db].append(hit_with_genome)
        for db, hits in db_results.items():
            if hits:
                summary_file = os.path.join(output_base, f"enteromark_{db}_abricate_summary.tsv")
                headers = list(hits[0].keys())
                with open(summary_file, 'w') as f:
                    f.write('\t'.join(headers) + '\n')
                    for hit in hits:
                        row = [str(hit.get(header, '')) for header in headers]
                        f.write('\t'.join(row) + '\n')
                self.logger.info("✓ Created %s summary: %s (%d hits)", db, summary_file, len(hits))
                self._create_database_summary_html(db, hits, output_base)
            else:
                self.logger.info("No hits for database %s, skipping summary", db)
    
    def _create_database_summary_html(self, database: str, hits: List[Dict], output_base: str):
        random_quote = self.get_random_quote()
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        unique_genomes = list(set(hit['genome'] for hit in hits))
        genes_per_genome = {}
        for hit in hits:
            genome = hit['genome']
            if genome not in genes_per_genome:
                genes_per_genome[genome] = set()
            genes_per_genome[genome].add(hit['gene'])
        
        html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroMark - Batch ABRicate Summary ({database.upper()})</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        body {{
            background: linear-gradient(135deg, #92400e 0%, #d97706 50%, #f59e0b 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
        }}
        .header {{
            text-align: center;
            margin-bottom: 30px;
        }}
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
            font-size: 10px;
            line-height: 1.1;
            white-space: pre;
            color: #fbbf24;
            text-shadow: 0 0 10px rgba(251, 191, 36, 0.5);
            overflow-x: auto;
        }}
        .quote-container {{
            background: rgba(255, 255, 255, 0.1);
            backdrop-filter: blur(10px);
            padding: 20px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            min-height: 100px;
            display: flex;
            flex-direction: column;
            justify-content: center;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: opacity 0.5s ease-in-out;
        }}
        .quote-text {{
            font-size: 18px;
            font-style: italic;
            margin-bottom: 10px;
        }}
        .quote-author {{
            font-size: 14px;
            color: #fbbf24;
            font-weight: bold;
        }}
        .report-section {{
            background: rgba(255, 255, 255, 0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 10px;
            margin-bottom: 20px;
            box-shadow: 0 4px 15px rgba(0, 0, 0, 0.2);
        }}
        .report-section h2 {{
            color: #b45309;
            border-bottom: 3px solid #d97706;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 24px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-top: 15px;
        }}
        .metric-card {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        .metric-label {{
            font-size: 14px;
            opacity: 0.9;
            margin-bottom: 5px;
        }}
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
        }}
        .table-responsive {{
            width: 100%;
            overflow-x: auto;
            margin: 20px 0;
        }}
        .summary-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 14px;
        }}
        .summary-table th {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }}
        .summary-table td {{
            padding: 12px;
            border-bottom: 1px solid #e5e7eb;
        }}
        .summary-table tr:nth-child(even) {{
            background-color: #fef3c7;
        }}
        .summary-table tr:hover {{
            background-color: #fffbeb;
        }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 20px;
            background: rgba(0, 0, 0, 0.3);
            border-radius: 10px;
            font-size: 14px;
        }}
        .timestamp {{
            color: #fbbf24;
            font-weight: bold;
        }}
        .authorship {{
            margin-top: 15px;
            padding: 15px;
            background: rgba(255, 255, 255, 0.1);
            border-radius: 8px;
            font-size: 12px;
        }}
        @media (max-width: 768px) {{
            .ascii-art {{ font-size: 6px; }}
            .metrics-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="ascii-container">
                <div class="ascii-art">

╔═══════════════════════════════════════════════════════════════════════════════╗
║   ███████╗███╗   ██╗████████╗███████╗██████╗  ██████╗ ███╗   ███╗ █████╗ ██████╗ ██╗  ██╗
║   ██╔════╝████╗  ██║╚══██╔══╝██╔════╝██╔══██╗██╔═══██╗████╗ ████║██╔══██╗██╔══██╗██║ ██╔╝
║   █████╗  ██╔██╗ ██║   ██║   █████╗  ██████╔╝██║   ██║██╔████╔██║███████║██████╔╝█████╔╝ 
║   ██╔══╝  ██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗██║   ██║██║╚██╔╝██║██╔══██║██╔══██╗██╔═██╗ 
║   ███████╗██║ ╚████║   ██║   ███████╗██║  ██║╚██████╔╝██║ ╚═╝ ██║██║  ██║██║  ██║██║  ██╗
║   ╚══════╝╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝
║                     EnteroMark - E. faecium ABRicate Analysis                 ║
╚═══════════════════════════════════════════════════════════════════════════════╝
                </div>
            </div>
            <div class="quote-container" id="quoteContainer">
                <div class="quote-text" id="quoteText">"{random_quote['text']}"</div>
                <div class="quote-author" id="quoteAuthor">— {random_quote['author']}</div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📊 Batch ABRicate Summary - {database.upper()} Database</h2>
            <div class="metrics-grid">
                <div class="metric-card">
                    <div class="metric-value">{len(hits)}</div>
                    <div class="metric-label">Total Hits</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{len(unique_genomes)}</div>
                    <div class="metric-label">Genomes</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{len(set(h['gene'] for h in hits))}</div>
                    <div class="metric-label">Unique Genes</div>
                </div>
                <div class="metric-card">
                    <div class="metric-value">{current_time.split()[0]}</div>
                    <div class="metric-label">Date</div>
                </div>
            </div>
        </div>
        
        <div class="report-section">
            <h2>🔍 Genes by Genome</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Genome</th>
                            <th>Count</th>
                            <th>Genes</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        for genome in sorted(unique_genomes):
            genes = genes_per_genome.get(genome, set())
            gene_list = ", ".join(sorted(genes))
            html_content += f"""
                    <tr>
                        <td><strong>{genome}</strong></td>
                        <td>{len(genes)}</td>
                        <td>{gene_list}</td>
                    </tr>
"""
        html_content += """
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="report-section">
            <h2>📈 Gene Frequency</h2>
            <div class="table-responsive">
                <table class="summary-table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Frequency</th>
                            <th>Genomes</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        gene_frequency = {}
        for hit in hits:
            gene = hit['gene']
            if gene not in gene_frequency:
                gene_frequency[gene] = set()
            gene_frequency[gene].add(hit['genome'])
        for gene, genomes in sorted(gene_frequency.items(), key=lambda x: len(x[1]), reverse=True):
            genome_list = ", ".join(sorted(genomes))
            html_content += f"""
                    <tr>
                        <td><strong>{gene}</strong></td>
                        <td>{len(genomes)}</td>
                        <td>{genome_list}</td>
                    </tr>
"""
        html_content += f"""
                    </tbody>
                </table>
            </div>
        </div>
        
        <div class="footer">
            <p><strong>ENTEROMARK</strong> - Batch ABRicate Summary</p>
            <p class="timestamp">Generated: {current_time}</p>
            <div class="authorship">
                <p><strong>Technical Support:</strong> Brown Beckley | GitHub: bbeckley-hub</p>
                <p>Email: brownbeckley94@gmail.com</p>
                <p>University of Ghana Medical School</p>
            </div>
        </div>
    </div>
    <script>
        const quotes = {json.dumps(self.science_quotes)};
        const quoteContainer = document.getElementById('quoteContainer');
        const quoteText = document.getElementById('quoteText');
        const quoteAuthor = document.getElementById('quoteAuthor');
        function getRandomQuote() {{ return quotes[Math.floor(Math.random() * quotes.length)]; }}
        function displayQuote() {{
            quoteContainer.style.opacity = '0';
            setTimeout(() => {{
                const quote = getRandomQuote();
                quoteText.textContent = '"' + quote.text + '"';
                quoteAuthor.textContent = '— ' + quote.author;
                quoteContainer.style.opacity = '1';
            }}, 500);
        }}
        setInterval(displayQuote, 10000);
    </script>
</body>
</html>"""
        html_file = os.path.join(output_base, f"enteromark_{database}_summary_report.html")
        with open(html_file, 'w') as f:
            f.write(html_content)
        self.logger.info("Database summary HTML report: %s", html_file)
    
    def create_database_json_summaries(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating JSON database summaries...")
        db_results = {}
        for genome_name, genome_result in all_results.items():
            for db, db_result in genome_result['results'].items():
                if db not in db_results:
                    db_results[db] = {
                        'hits': [],
                        'genomes': [],
                        'gene_frequency': {}
                    }
                for hit in db_result['hits']:
                    hit_with_genome = hit.copy()
                    hit_with_genome['genome'] = genome_name
                    db_results[db]['hits'].append(hit_with_genome)
                if genome_name not in db_results[db]['genomes']:
                    db_results[db]['genomes'].append(genome_name)
        for db, data in db_results.items():
            if not data['hits']:
                continue
            gene_frequency = {}
            for hit in data['hits']:
                gene = hit['gene']
                if gene not in gene_frequency:
                    gene_frequency[gene] = {
                        'count': 0,
                        'genomes': set(),
                        'details': []
                    }
                gene_frequency[gene]['count'] += 1
                gene_frequency[gene]['genomes'].add(hit['genome'])
                gene_frequency[gene]['details'].append({
                    'genome': hit['genome'],
                    'product': hit['product'],
                    'coverage': hit['coverage_percent'],
                    'identity': hit['identity_percent'],
                    'accession': hit['accession']
                })
            for gene in gene_frequency:
                gene_frequency[gene]['genomes'] = list(gene_frequency[gene]['genomes'])
            json_summary = {
                'metadata': {
                    'database': db,
                    'analysis_date': self.metadata['analysis_date'],
                    'tool': self.metadata['tool_name'],
                    'version': self.metadata['version'],
                    'total_hits': len(data['hits']),
                    'total_genomes': len(data['genomes']),
                    'unique_genes': len(gene_frequency)
                },
                'gene_frequency': gene_frequency,
                'hits': data['hits'][:100000]
            }
            json_file = os.path.join(output_base, f"enteromark_{db}_summary.json")
            with open(json_file, 'w') as f:
                json.dump(json_summary, f, indent=2, default=str)
            self.logger.info("✓ Created JSON summary: %s", json_file)
    
    def create_master_json_summary(self, all_results: Dict[str, Any], output_base: str):
        self.logger.info("Creating master JSON summary...")
        master_summary = {
            'metadata': {
                'tool': 'EnteroMark ABRicate Module',
                'version': self.metadata['version'],
                'analysis_date': self.metadata['analysis_date'],
                'total_genomes': len(all_results),
                'databases_used': self.required_databases
            },
            'genome_summaries': {},
            'critical_findings': {},
            'cross_database_patterns': {}
        }
        all_hits_by_gene = {}
        for genome_name, genome_result in all_results.items():
            all_genome_hits = []
            for db_result in genome_result['results'].values():
                all_genome_hits.extend(db_result['hits'])
            analysis = self.analyze_efaecalis_genes(all_genome_hits)
            master_summary['genome_summaries'][genome_name] = {
                'total_hits': genome_result['total_hits'],
                'critical_resistance': len(analysis['critical_resistance_genes']),
                'high_risk_virulence': len(analysis['high_risk_virulence_genes']),
                'beta_lactamases': len(analysis['beta_lactamase_genes']),
                'other_genes': len(analysis['other_genes'])
            }
            if analysis['critical_resistance_genes']:
                if 'critical_genomes' not in master_summary['critical_findings']:
                    master_summary['critical_findings']['critical_genomes'] = []
                master_summary['critical_findings']['critical_genomes'].append(genome_name)
            if analysis['high_risk_virulence_genes']:
                if 'hv_genomes' not in master_summary['critical_findings']:
                    master_summary['critical_findings']['hv_genomes'] = []
                master_summary['critical_findings']['hv_genomes'].append(genome_name)
            for hit in all_genome_hits:
                gene = hit['gene']
                if gene not in all_hits_by_gene:
                    all_hits_by_gene[gene] = {
                        'count': 0,
                        'genomes': set(),
                        'products': set(),
                        'databases': set()
                    }
                all_hits_by_gene[gene]['count'] += 1
                all_hits_by_gene[gene]['genomes'].add(genome_name)
                all_hits_by_gene[gene]['products'].add(hit['product'])
                all_hits_by_gene[gene]['databases'].add(hit['database'])
        for gene in all_hits_by_gene:
            all_hits_by_gene[gene]['genomes'] = list(all_hits_by_gene[gene]['genomes'])
            all_hits_by_gene[gene]['products'] = list(all_hits_by_gene[gene]['products'])
            all_hits_by_gene[gene]['databases'] = list(all_hits_by_gene[gene]['databases'])
        master_summary['cross_database_patterns'] = {
            'total_genes_found': len(all_hits_by_gene),
            'common_genes': {g: d for g, d in all_hits_by_gene.items() if d['count'] > 1},
            'top_genes': sorted(
                [(g, d) for g, d in all_hits_by_gene.items()],
                key=lambda x: x[1]['count'],
                reverse=True
            )[:50]
        }
        json_file = os.path.join(output_base, "enteromark_abricate_master_summary.json")
        with open(json_file, 'w') as f:
            json.dump(master_summary, f, indent=2, default=str)
        self.logger.info("✓ Created master JSON summary: %s", json_file)
    
    def process_single_genome(self, genome_file: str, output_base: str = "abricate_results") -> Dict[str, Any]:
        genome_name = Path(genome_file).stem
        results_dir = os.path.join(output_base, genome_name)
        os.makedirs(results_dir, exist_ok=True)
        databases = self.required_databases
        results = {}
        for db in databases:
            result = self.run_abricate_single_db(genome_file, db, results_dir)
            results[db] = result
        self.create_comprehensive_html_report(genome_name, results, results_dir)
        return {
            'genome': genome_name,
            'results': results,
            'total_hits': sum(r['hit_count'] for r in results.values())
        }
    
    def process_multiple_genomes(self, genome_pattern: str, output_base: str = "abricate_results") -> Dict[str, Any]:
        if not self.check_abricate_installed():
            raise RuntimeError("ABRicate not properly installed")
        self.setup_abricate_databases()
        fasta_patterns = [genome_pattern, f"{genome_pattern}.fasta", f"{genome_pattern}.fa",
                          f"{genome_pattern}.fna", f"{genome_pattern}.faa"]
        genome_files = []
        for pattern in fasta_patterns:
            genome_files.extend(glob.glob(pattern))
        genome_files = list(set(genome_files))
        if not genome_files:
            raise FileNotFoundError(f"No FASTA files found matching pattern: {genome_pattern}")
        self.logger.info("Found %d genomes", len(genome_files))
        os.makedirs(output_base, exist_ok=True)
        all_results = {}
        if len(genome_files) > 1 and self.cpus > 1:
            with ThreadPoolExecutor(max_workers=self.cpus) as executor:
                future_to_genome = {executor.submit(self.process_single_genome, g, output_base): g for g in genome_files}
                for future in as_completed(future_to_genome):
                    genome = future_to_genome[future]
                    try:
                        result = future.result()
                        all_results[Path(genome).stem] = result
                        self.logger.info("✓ Completed: %s (%d hits)", result['genome'], result['total_hits'])
                    except Exception as e:
                        self.logger.error("✗ Failed: %s - %s", genome, e)
        else:
            for genome in genome_files:
                try:
                    result = self.process_single_genome(genome, output_base)
                    all_results[Path(genome).stem] = result
                    self.logger.info("✓ Completed: %s (%d hits)", result['genome'], result['total_hits'])
                except Exception as e:
                    self.logger.error("✗ Failed: %s - %s", genome, e)
        self.create_database_summaries(all_results, output_base)
        self.create_database_json_summaries(all_results, output_base)
        self.create_master_json_summary(all_results, output_base)
        self.logger.info("=== ANALYSIS COMPLETE ===")
        return all_results

def main():
    parser = argparse.ArgumentParser(
        description='EnteroMark ABRicate Analysis - MAXIMUM SPEED VERSION',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  python enteromark_abricate.py "*.fasta"
  python enteromark_abricate.py "*.fna" --output my_results --cpus 8"""
    )
    parser.add_argument('pattern', help='File pattern for genomes (e.g., "*.fasta")')
    parser.add_argument('--cpus', '-c', type=int, default=None, help='Number of CPU cores (default: auto-detect)')
    parser.add_argument('--output', '-o', default='abricate_results', help='Output directory')
    args = parser.parse_args()
    executor = EnteroMarkAbricateExecutor(cpus=args.cpus)
    try:
        results = executor.process_multiple_genomes(args.pattern, args.output)
        executor.logger.info("\n" + "="*50)
        executor.logger.info("📊 FINAL SUMMARY")
        executor.logger.info("="*50)
        total_crit = 0
        total_high = 0
        for genome_name, result in results.items():
            all_hits = []
            for db_result in result['results'].values():
                all_hits.extend(db_result['hits'])
            analysis = executor.analyze_efaecalis_genes(all_hits)
            executor.logger.info("✓ %s: %d hits (critical: %d, high-risk virulence: %d, beta-lactamases: %d)",
                               genome_name, result['total_hits'],
                               analysis['total_critical_resistance'],
                               analysis['total_high_risk_virulence'],
                               analysis['total_beta_lactamase'])
            total_crit += analysis['total_critical_resistance']
            total_high += analysis['total_high_risk_virulence']
        executor.logger.info("\n📁 Results saved to: %s", args.output)
        executor.logger.info("⚡ MAXIMUM SPEED: %d cores utilized", executor.cpus)
    except Exception as e:
        executor.logger.error("Analysis failed: %s", e)
        sys.exit(1)

if __name__ == "__main__":
    main()