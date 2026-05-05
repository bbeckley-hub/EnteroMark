#!/usr/bin/env python3
"""
EnteroMark Ultimate Summary Module – Gene‑Centric, Interactive HTML Report for E. faecium
Enhanced with Key Statistics, AI Guide, Call to Action, Detailed Info, and Improved UX
Author: Brown Beckley <brownbeckley94@gmail.com>
Affiliation: University of Ghana Medical School
Version: 1.0.0 – Fixed search bars and QC header sorting
Date: 2026-05-05
"""

import os
import sys
import json
import re
import argparse
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any
from datetime import datetime
from collections import defaultdict, Counter
import warnings
warnings.filterwarnings('ignore')

from bs4 import BeautifulSoup
import re

# =============================================================================
# HTML PARSER – parses all EnteroMark summary reports
# =============================================================================
class EnteroMarkHTMLParser:
    def __init__(self):
        self.db_files = [
            'enteromark_card_summary_report.html',
            'enteromark_resfinder_summary_report.html',
            'enteromark_argannot_summary_report.html',
            'enteromark_vfdb_summary_report.html',
            'enteromark_plasmidfinder_summary_report.html',
            'enteromark_megares_summary_report.html',
            'enteromark_ncbi_summary_report.html',
            'enteromark_bacmet2_summary_report.html'
        ]

    def normalize_sample_id(self, sample_id: str) -> str:
        sample = str(sample_id).strip()
        for ext in ['.fna', '.fasta', '.fa', '.gb', '.txt', '.tsv', '.fna.gz']:
            if sample.endswith(ext):
                sample = sample[:-len(ext)]
        return Path(sample).name

    def parse_qc_summary(self, file_path: Path) -> Dict[str, Dict]:
        """Parse enteromark_fasta_qc_summary.html by extracting embedded JSON data."""
        print(f"  🧬 Parsing FASTA QC: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            match = re.search(r'const\s+tableData\s*=\s*(\[[\s\S]*?\]);', html)
            if not match:
                match = re.search(r'var\s+tableData\s*=\s*(\[[\s\S]*?\]);', html)
            if not match:
                print("    ⚠️ Could not find tableData JSON")
                return {}
            data_list = json.loads(match.group(1))
            results = {}
            for item in data_list:
                sample_raw = item.get('filename', '')
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                qc_data = {}
                key_map = {
                    'total_sequences': 'total_sequences',
                    'total_length': 'total_length',
                    'total_bases': 'total_bases',
                    'gc_percent': 'gc_percent',
                    'at_percent': 'at_percent',
                    'n50': 'n50',
                    'n75': 'n75',
                    'n90': 'n90',
                    'median_length': 'median_length',
                    'mean_length': 'mean_length',
                    'longest_sequence': 'longest_sequence',
                    'shortest_sequence': 'shortest_sequence',
                    'ambiguous_percent': 'ambiguous_percent',
                    'sequences_with_n': 'sequences_with_n',
                    'max_n_run': 'max_n_run',
                    'homopolymers_count': 'homopolymers_count',
                    'max_homopolymer': 'max_homopolymer',
                    'duplicate_sequences': 'duplicate_sequences',
                    'short_sequences': 'short_sequences',
                    'long_sequences': 'long_sequences',
                    'file_size_mb': 'file_size_mb',
                    'warnings_count': 'warnings',
                    'status_str': 'ef_status'
                }
                for json_key, internal_key in key_map.items():
                    val = item.get(json_key)
                    if val is None or val == '':
                        qc_data[internal_key] = 'ND'
                    else:
                        qc_data[internal_key] = val
                results[sample] = qc_data
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error: {e}")
            return {}

    def parse_mlst_summary(self, file_path: Path) -> Dict[str, Dict]:
        print(f"  🧬 Parsing MLST: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            table = soup.find('table', id='mlst-summary-table')
            if not table:
                table = soup.find('table', class_='summary-table')
            if not table:
                return {}
            rows = table.find_all('tr')
            if len(rows) < 2:
                return {}
            results = {}
            for row in rows[1:]:
                cols = row.find_all('td')
                if len(cols) < 2:
                    continue
                sample_raw = cols[0].get_text().strip()
                if not sample_raw:
                    continue
                sample = self.normalize_sample_id(sample_raw)
                st_raw = cols[1].get_text().strip().replace('ST', '')
                st = st_raw if st_raw and st_raw.lower() != 'nd' else 'ND'
                results[sample] = {'ST': st}
            print(f"    ✓ Parsed {len(results)} samples")
            return results
        except Exception as e:
            print(f"    ❌ Error: {e}")
            return {}

    def parse_abricate_database_summary(self, file_path: Path, total_samples: int = 0) -> Dict[str, Dict]:
        print(f"  🧬 Parsing DB: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}
            db_name = file_path.stem.replace('enteromark_', '').replace('_summary_report', '')
            if db_name.endswith('_summary'):
                db_name = db_name.replace('_summary', '')
            df_freq = self._parse_html_table_to_dict(tables[1])
            gene_freq = {}
            for row in df_freq:
                gene_full = row.get('Gene', '').strip()
                if not gene_full:
                    continue
                gene = re.sub(r'^\([^)]+\)', '', gene_full).strip()
                if not gene:
                    gene = gene_full
                freq_str = row.get('Frequency', '0')
                count_match = re.search(r'(\d+)', freq_str)
                count = int(count_match.group(1)) if count_match else 0
                percentage = (count / total_samples * 100) if total_samples else 0
                genomes_str = row.get('Genomes', '')
                genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()] if genomes_str else []
                gene_freq[gene] = {
                    'count': count,
                    'percentage': round(percentage, 2),
                    'frequency_display': f"{count} ({percentage:.1f}%)",
                    'genomes': genomes,
                    'database': db_name
                }
            print(f"    ✓ {db_name.upper()}: {len(gene_freq)} genes")
            return gene_freq
        except Exception as e:
            print(f"    ❌ Error: {e}")
            return {}

    def parse_amrfinder_summary(self, file_path: Path, total_samples: int = 0) -> Dict[str, Dict]:
        print(f"  🧬 Parsing AMRfinder: {file_path.name}")
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                html = f.read()
            soup = BeautifulSoup(html, 'html.parser')
            tables = soup.find_all('table')
            if len(tables) < 2:
                return {}
            df_freq = self._parse_html_table_to_dict(tables[1])
            gene_freq = {}
            for row in df_freq:
                gene = row.get('Gene', '').strip()
                if not gene:
                    continue
                freq_str = row.get('Frequency', '0')
                count_match = re.search(r'(\d+)', freq_str)
                count = int(count_match.group(1)) if count_match else 0
                percentage = (count / total_samples * 100) if total_samples else 0
                genomes_str = row.get('Genomes', '')
                genomes = [self.normalize_sample_id(g.strip()) for g in genomes_str.split(',') if g.strip()] if genomes_str else []
                gene_freq[gene] = {
                    'count': count,
                    'percentage': round(percentage, 2),
                    'frequency_display': f"{count} ({percentage:.1f}%)",
                    'genomes': genomes,
                    'database': 'amrfinder'
                }
            print(f"    ✓ AMRFinder: {len(gene_freq)} genes")
            return gene_freq
        except Exception as e:
            print(f"    ❌ Error: {e}")
            return {}

    def _parse_html_table_to_dict(self, table) -> List[Dict]:
        rows = table.find_all('tr')
        if not rows:
            return []
        headers = [th.get_text().strip() for th in rows[0].find_all(['th', 'td'])]
        data = []
        for row in rows[1:]:
            cols = row.find_all(['td', 'th'])
            row_dict = {}
            for i, col in enumerate(cols):
                if i < len(headers):
                    row_dict[headers[i]] = col.get_text().strip()
            data.append(row_dict)
        return data


# =============================================================================
# DATA ANALYZER (enhanced for E. faecium)
# =============================================================================
class EnteroMarkDataAnalyzer:
    def __init__(self):
        self.critical_resistance = {
            'vanA', 'vanB', 'vanD', 'vanM', 'optrA', 'cfr', 'cfrB', 'poxtA',
            'aac(6\')-Ii', 'ant(6)-Ia', 'aph(3\')-III', 'aac(6\')-Ie-aph(2\'\')-Ia'
        }
        self.high_risk_virulence = {
            'esp', 'ace', 'asal', 'gelE', 'cylA', 'cylB', 'cylL', 'cylM', 'cylR',
            'efaA', 'fsrA', 'fsrB', 'fsrC', 'hyl', 'scm', 'acm', 'ecbA', 'ecbB'
        }
        self.beta_lactamase = {'pbp5', 'pbp5_S462A', 'pbp5_E629V', 'pbp5_M485A', 'pbp5_M426I'}

    def categorize_gene(self, gene: str) -> str:
        g = gene.lower()
        if any(crit.lower() in g for crit in self.critical_resistance):
            return 'Critical Resistance'
        elif any(vir.lower() in g for vir in self.high_risk_virulence):
            return 'High-Risk Virulence'
        elif any(bla.lower() in g for bla in self.beta_lactamase):
            return 'Beta-Lactamase'
        else:
            return 'Other'

    def create_patterns(self, samples_data: Dict, gene_freqs: Dict) -> Dict:
        patterns = {
            'st_distribution': Counter(),
            'st_samples': defaultdict(list),
            'high_risk_combinations': [],
            'van_plus_pbp5': [],
            'linezolid_plus_van': []
        }
        sample_genes = defaultdict(set)
        for db, genes in gene_freqs.items():
            for gene, data in genes.items():
                for genome in data['genomes']:
                    sample_genes[genome].add(gene)
        for sample, data in samples_data.items():
            st = data.get('mlst', {}).get('ST', 'ND')
            if st != 'ND':
                patterns['st_distribution'][st] += 1
                patterns['st_samples'][st].append(sample)
            genes = sample_genes.get(sample, set())
            crit = [g for g in genes if self.categorize_gene(g) == 'Critical Resistance']
            hv = [g for g in genes if self.categorize_gene(g) == 'High-Risk Virulence']
            if crit and hv:
                patterns['high_risk_combinations'].append({
                    'sample': sample,
                    'st': st,
                    'critical_resistance': crit,
                    'high_risk_virulence': hv
                })
            van_genes = [g for g in genes if any(v in g.lower() for v in ['vana', 'vanb', 'vand', 'vanm'])]
            pbp5_genes = [g for g in genes if 'pbp5' in g.lower()]
            if van_genes and pbp5_genes:
                patterns['van_plus_pbp5'].append({
                    'sample': sample,
                    'st': st,
                    'van_genes': van_genes,
                    'pbp5_genes': pbp5_genes
                })
            linezolid_genes = [g for g in genes if any(v in g.lower() for v in ['optr', 'cfr', 'poxt'])]
            if linezolid_genes and van_genes:
                patterns['linezolid_plus_van'].append({
                    'sample': sample,
                    'st': st,
                    'linezolid_genes': linezolid_genes,
                    'van_genes': van_genes
                })
        return patterns


# =============================================================================
# ENHANCED HTML GENERATOR
# =============================================================================
class EnteroMarkHTMLGenerator:
    def __init__(self, analyzer: EnteroMarkDataAnalyzer):
        self.analyzer = analyzer

    def generate_main_report(self, integrated_data: Dict, output_dir: Path) -> str:
        print("\n🎨 Generating EnteroMark Ultimate HTML report...")
        html = self._create_ultimate_html(integrated_data)
        output_file = output_dir / "enteromark_ultimate_report.html"
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        print(f"    ✅ HTML report saved: {output_file}")
        return str(output_file)

    def _create_ultimate_html(self, data: Dict) -> str:
        metadata = data.get('metadata', {})
        samples = data.get('samples', {})
        qc_data = data.get('qc_data', {})
        mlst_data = data.get('mlst_data', {})
        gene_freqs = data.get('gene_frequencies', {})
        patterns = data.get('patterns', {})
        total_samples = len(samples)

        # Compute dashboard numbers
        unique_sts = len(patterns.get('st_distribution', {}))
        amr_genes = sum(1 for db, genes in gene_freqs.items() if db not in ['vfdb', 'plasmidfinder'] for _ in genes)
        virulence_genes = len(gene_freqs.get('vfdb', {}))
        plasmid_genes = len(gene_freqs.get('plasmidfinder', {}))
        high_risk_combos = len(patterns.get('high_risk_combinations', []))

        # Sample overview rows
        sample_rows = []
        for sample in sorted(samples.keys()):
            qc = qc_data.get(sample, {})
            mlst = mlst_data.get(sample, {})
            st = mlst.get('ST', 'ND')
            n50 = qc.get('n50', 'ND')
            gc = qc.get('gc_percent', 'ND')
            at = qc.get('at_percent', 'ND')
            n75 = qc.get('n75', 'ND')
            n90 = qc.get('n90', 'ND')
            total_bases = qc.get('total_bases', 'ND')
            sample_rows.append(f"""
            <tr>
                <td class="col-sample">{sample}</td>
                <td class="col-st">{st}</td>
                <td class="col-n50">{n50}</td>
                <td class="col-gc">{gc}</td>
                <td class="col-at">{at}</td>
                <td class="col-n75">{n75}</td>
                <td class="col-n90">{n90}</td>
                <td class="col-total-bases">{total_bases}</td>
            </tr>
            """)

        # AMR table
        amr_rows = []
        for db, genes in gene_freqs.items():
            if db in ['vfdb', 'plasmidfinder']:
                continue
            for gene, gdata in genes.items():
                genome_tags = ''.join([f'<span class="genome-tag">{gen}</span>' for gen in gdata['genomes']])
                amr_rows.append(f"""
                <tr>
                    <td class="col-gene"><strong>{gene}</strong></td>
                    <td class="col-database">{db.upper()}</td>
                    <td class="col-frequency">{gdata['frequency_display']}</td>
                    <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                </tr>
                """)
        amr_rows.sort(key=lambda x: int(re.search(r'(\d+)', x.split('<td class="col-frequency">')[1].split('<')[0])[0] if re.search(r'(\d+)', x) else 0), reverse=True)

        # Virulence table
        virulence_rows = []
        if 'vfdb' in gene_freqs:
            for gene, gdata in sorted(gene_freqs['vfdb'].items(), key=lambda x: x[1]['count'], reverse=True):
                genome_tags = ''.join([f'<span class="genome-tag">{gen}</span>' for gen in gdata['genomes']])
                virulence_rows.append(f"""
                <tr>
                    <td class="col-gene"><strong>{gene}</strong></td>
                    <td class="col-database">VFDB</td>
                    <td class="col-frequency">{gdata['frequency_display']}</td>
                    <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                </tr>
                """)

        # Plasmid table
        plasmid_rows = []
        if 'plasmidfinder' in gene_freqs:
            for gene, gdata in sorted(gene_freqs['plasmidfinder'].items(), key=lambda x: x[1]['count'], reverse=True):
                genome_tags = ''.join([f'<span class="genome-tag">{gen}</span>' for gen in gdata['genomes']])
                plasmid_rows.append(f"""
                <tr>
                    <td class="col-gene"><strong>{gene}</strong></td>
                    <td class="col-database">PlasmidFinder</td>
                    <td class="col-frequency">{gdata['frequency_display']}</td>
                    <td class="col-genomes"><div class="genome-list">{genome_tags}</div></td>
                </tr>
                """)

        # QC metrics headers
        all_metrics = set()
        for qc in qc_data.values():
            all_metrics.update(qc.keys())
        col_order = ['total_sequences', 'total_length', 'total_bases', 'gc_percent', 'at_percent',
                     'n50', 'n75', 'n90', 'median_length', 'mean_length', 'longest_sequence',
                     'shortest_sequence', 'ambiguous_percent', 'sequences_with_n', 'max_n_run',
                     'homopolymers_count', 'max_homopolymer', 'duplicate_sequences',
                     'short_sequences', 'long_sequences', 'file_size_mb', 'warnings', 'ef_status']
        qc_headers = [h for h in col_order if h in all_metrics]
        qc_rows = []
        for sample, qc in sorted(qc_data.items()):
            row = f"<tr><td class='col-sample'>{sample}"
            for metric in qc_headers:
                val = qc.get(metric, 'ND')
                row += f"<td class='col-{metric}'>{val}"
            row += "</tr>"
            qc_rows.append(row)

        # MLST distribution rows with samples
        st_rows = ""
        total_st = sum(patterns.get('st_distribution', {}).values())
        for st, cnt in sorted(patterns.get('st_distribution', {}).items(), key=lambda x: x[1], reverse=True):
            pct = cnt / total_st * 100 if total_st else 0
            samples_list = ', '.join(patterns.get('st_samples', {}).get(st, []))
            st_rows += f"<tr><td class='col-st'>ST{st}</td><td class='col-count'>{cnt}</td><td class='col-pct'>{pct:.1f}%</td><td class='col-samples'><div class='genome-list' style='max-height:100px; overflow-y:auto;'>{samples_list}</div></td></tr>"

        # Pattern tables
        high_risk_rows = ""
        for combo in patterns.get('high_risk_combinations', [])[:20]:
            high_risk_rows += f"<tr><td>{combo['sample']}</td><td>{combo['st']}</td><td>{', '.join(combo['critical_resistance'])}</td><td>{', '.join(combo['high_risk_virulence'])}</td></tr>"
        van_pbp5_rows = ""
        for combo in patterns.get('van_plus_pbp5', [])[:20]:
            van_pbp5_rows += f"<tr><td>{combo['sample']}</td><td>{combo['st']}</td><td>{', '.join(combo['van_genes'])}</td><td>{', '.join(combo['pbp5_genes'])}</td></tr>"
        linezolid_van_rows = ""
        for combo in patterns.get('linezolid_plus_van', [])[:20]:
            linezolid_van_rows += f"<tr><td>{combo['sample']}</td><td>{combo['st']}</td><td>{', '.join(combo['linezolid_genes'])}</td><td>{', '.join(combo['van_genes'])}</td></tr>"

        # JavaScript with fixed search and sort functions
        js = f"""
        <script>
            function filterTable(tableId, inputId) {{
                const filter = document.getElementById(inputId).value.toUpperCase();
                const table = document.getElementById(tableId);
                const rows = table.getElementsByTagName('tr');
                for (let i = 1; i < rows.length; i++) {{
                    let txt = rows[i].textContent.toUpperCase();
                    rows[i].style.display = txt.indexOf(filter) > -1 ? '' : 'none';
                }}
            }}

            function highlightGenomes(tableId, inputId) {{
                const filter = document.getElementById(inputId).value.trim().toUpperCase();
                const table = document.getElementById(tableId);
                const allTags = table.querySelectorAll('.genome-tag');
                allTags.forEach(tag => tag.classList.remove('highlight-genome'));
                if (filter === "") return;
                allTags.forEach(tag => {{
                    if (tag.textContent.toUpperCase().includes(filter)) {{
                        tag.classList.add('highlight-genome');
                    }}
                }});
            }}

            function exportTableToCSV(tableId, filename) {{
                const table = document.getElementById(tableId);
                let csv = [];
                const rows = table.querySelectorAll('tr');
                for (let row of rows) {{
                    const cols = row.querySelectorAll('td,th');
                    const rowData = Array.from(cols).map(c => '"' + c.innerText.replace(/"/g, '""') + '"');
                    csv.push(rowData.join(','));
                }}
                const blob = new Blob([csv.join('\\n')], {{type: 'text/csv'}});
                const a = document.createElement('a');
                a.href = URL.createObjectURL(blob);
                a.download = filename;
                a.click();
                URL.revokeObjectURL(a.href);
            }}

            function sortTable(tableId, colIndex, type) {{
                const table = document.getElementById(tableId);
                const tbody = table.tBodies[0];
                const rows = Array.from(tbody.rows);
                const isAscending = table.getAttribute('data-sort-dir') !== 'asc';
                rows.sort((a, b) => {{
                    let aVal = a.cells[colIndex].innerText.trim();
                    let bVal = b.cells[colIndex].innerText.trim();
                    if (type === 'number') {{
                        aVal = parseFloat(aVal.replace(/,/g, '')) || 0;
                        bVal = parseFloat(bVal.replace(/,/g, '')) || 0;
                        return isAscending ? aVal - bVal : bVal - aVal;
                    }} else {{
                        return isAscending ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                    }}
                }});
                tbody.append(...rows);
                table.setAttribute('data-sort-dir', isAscending ? 'asc' : 'desc');
            }}

            function addSortListeners(tableId, colTypes) {{
                const table = document.getElementById(tableId);
                if (!table) return;
                const headers = table.querySelectorAll('th');
                headers.forEach((header, idx) => {{
                    if (colTypes[idx]) {{
                        header.style.cursor = 'pointer';
                        header.addEventListener('click', () => sortTable(tableId, idx, colTypes[idx]));
                    }}
                }});
            }}

            function setFilter(term, tableId, searchInputId) {{
                document.getElementById(searchInputId).value = term;
                filterTable(tableId, searchInputId);
            }}

            function switchTab(tabName) {{
                document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
                document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
                document.getElementById(tabName + '-tab').classList.add('active');
                event.currentTarget.classList.add('active');
            }}

            // Initialize sort listeners after DOM loads
            document.addEventListener('DOMContentLoaded', function() {{
                // Samples table: first col string, rest numbers (8 columns total)
                addSortListeners('samples-table', ['string', 'string', 'number', 'number', 'number', 'number', 'number', 'number']);
                // MLST table: string, number, number, string
                addSortListeners('mlst-table', ['string', 'number', 'number', 'string']);
                // AMR table: string, string, number, string
                addSortListeners('amr-table', ['string', 'string', 'number', 'string']);
                // Virulence table: same as AMR
                addSortListeners('virulence-table', ['string', 'string', 'number', 'string']);
                // Plasmid table: same
                addSortListeners('plasmid-table', ['string', 'string', 'number', 'string']);
                // QC table: first column string, all others number
                const qcHeaders = document.querySelectorAll('#qc-table th');
                const qcColTypes = ['string'] + Array(qcHeaders.length - 1).fill('number');
                addSortListeners('qc-table', qcColTypes);
                // Show first tab
                document.querySelector('.tab-button').click();
            }});
        </script>
        """

        # Build HTML
        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>EnteroMark Ultimate Report – E. faecium</title>
    <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css">
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        body {{
            background: linear-gradient(135deg, #92400e 0%, #d97706 50%, #f59e0b 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #ffffff;
            padding: 20px;
            min-height: 100vh;
        }}
        .container {{ max-width: 1600px; margin: 0 auto; }}
        .main-header {{
            background: rgba(0,0,0,0.7);
            padding: 30px;
            border-radius: 15px;
            margin-bottom: 30px;
            text-align: center;
            border: 2px solid rgba(245,158,11,0.5);
        }}
        .main-header h1 {{ font-size: 2.8em; margin-bottom: 10px; color: #fbbf24; }}
        .metadata-bar {{
            background: rgba(255,255,255,0.1);
            padding: 15px;
            border-radius: 10px;
            margin: 20px 0;
            display: flex;
            justify-content: center;
            gap: 30px;
            flex-wrap: wrap;
        }}
        .dashboard-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        .dashboard-card {{
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 20px;
            border-radius: 12px;
            text-align: center;
            cursor: pointer;
            transition: transform 0.2s, box-shadow 0.2s;
            box-shadow: 0 4px 12px rgba(0,0,0,0.2);
            border-left: 5px solid #d97706;
        }}
        .dashboard-card:hover {{ transform: translateY(-5px); box-shadow: 0 8px 20px rgba(0,0,0,0.3); }}
        .card-number {{ font-size: 2.5em; font-weight: bold; color: #d97706; }}
        .tab-navigation {{
            display: flex;
            gap: 5px;
            flex-wrap: wrap;
            margin-bottom: 20px;
            background: rgba(0,0,0,0.5);
            padding: 15px;
            border-radius: 12px;
            position: sticky;
            top: 10px;
            z-index: 100;
        }}
        .tab-button {{
            padding: 10px 20px;
            background: #f5f5f5;
            border: none;
            border-radius: 8px;
            cursor: pointer;
            font-weight: bold;
            color: #333;
            transition: all 0.3s;
        }}
        .tab-button.active {{
            background: #d97706;
            color: white;
        }}
        .tab-content {{
            display: none;
            background: rgba(255,255,255,0.95);
            color: #1f2937;
            padding: 25px;
            border-radius: 12px;
            margin-bottom: 20px;
            animation: fadeIn 0.3s;
        }}
        .tab-content.active {{ display: block; }}
        @keyframes fadeIn {{ from {{ opacity: 0; }} to {{ opacity: 1; }} }}
        .section-header {{
            color: #b45309;
            border-bottom: 3px solid #d97706;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 1.8em;
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
        }}
        .info-box {{
            background: #e0f2fe;
            border-left: 5px solid #0284c7;
            padding: 15px;
            margin: 15px 0;
            border-radius: 8px;
            color: #0c4a6e;
        }}
        .info-box h4 {{ margin-bottom: 10px; color: #0369a1; }}
        .info-box ul {{ margin-left: 20px; }}
        .search-box {{
            width: 100%;
            padding: 10px;
            margin: 15px 0;
            border: 2px solid #ddd;
            border-radius: 8px;
            font-size: 1em;
        }}
        .action-buttons {{
            display: flex;
            flex-wrap: wrap;
            gap: 10px;
            margin: 15px 0;
        }}
        .action-btn {{
            padding: 8px 16px;
            background: #d97706;
            color: white;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-weight: bold;
            transition: background 0.2s;
        }}
        .action-btn:hover {{ background: #b45309; }}
        .master-scrollable-container {{
            overflow-x: auto;
            max-height: 600px;
            overflow-y: auto;
            border: 1px solid #ddd;
            border-radius: 8px;
        }}
        .data-table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 0.9em;
            min-width: 800px;
        }}
        .data-table th {{
            background: linear-gradient(135deg, #d97706 0%, #b45309 100%);
            color: white;
            padding: 12px;
            text-align: left;
            position: sticky;
            top: 0;
            cursor: pointer;
            user-select: none;
        }}
        .data-table td {{
            padding: 10px;
            border-bottom: 1px solid #e5e7eb;
            vertical-align: top;
        }}
        .data-table tr:nth-child(even) {{ background-color: #fef3c7; }}
        .data-table tr:hover {{ background-color: #fffbeb; }}
        .highlight-genome {{
            background-color: #a3e635 !important;
            font-weight: bold;
        }}
        .genome-list {{
            max-height: 150px;
            overflow-y: auto;
            display: flex;
            flex-wrap: wrap;
            gap: 4px;
        }}
        .genome-tag {{
            display: inline-block;
            background: #e0f2f1;
            color: #2c7a4d;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 0.75em;
            margin: 2px;
            border: 1px solid #b2dfdb;
            white-space: nowrap;
        }}
        .alert-box {{
            padding: 15px;
            border-radius: 8px;
            margin: 15px 0;
            display: flex;
            align-items: center;
            gap: 15px;
        }}
        .alert-info {{ background: #d1ecf1; color: #0c5460; border-left: 5px solid #17a2b8; }}
        .footer {{
            text-align: center;
            margin-top: 30px;
            padding: 30px;
            background: rgba(0,0,0,0.8);
            border-radius: 15px;
        }}
        .footer a {{ color: #ffc107; text-decoration: none; }}
        .footer a:hover {{ text-decoration: underline; }}
        .critical-cards {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        .critical-card {{
            background: #fff;
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-left: 4px solid;
        }}
        @media (max-width: 768px) {{
            .tab-button {{ padding: 6px 12px; font-size: 0.8em; }}
        }}
    </style>
    {js}
</head>
<body>
<div class="container">
    <div class="main-header">
        <h1><i class="fas fa-dna"></i> EnteroMark Ultimate Report</h1>
        <p>Comprehensive AMR, Virulence, MLST, and QC analysis for <strong>Enterococcus faecium</strong></p>
        <div class="metadata-bar">
            <span><i class="fas fa-calendar"></i> {metadata.get('analysis_date', 'Unknown')}</span>
            <span><i class="fas fa-database"></i> {total_samples} samples</span>
            <span><i class="fas fa-tools"></i> EnteroMark v1.0.0</span>
        </div>
    </div>

    <!-- Dashboard Cards -->
    <div class="dashboard-grid">
        <div class="dashboard-card" onclick="switchTab('samples')"><div class="card-number">{total_samples}</div><div>Samples</div><i class="fas fa-vial fa-2x" style="color:#d97706"></i></div>
        <div class="dashboard-card" onclick="switchTab('mlst')"><div class="card-number">{unique_sts}</div><div>Unique STs</div><i class="fas fa-code-branch fa-2x"></i></div>
        <div class="dashboard-card" onclick="switchTab('amr')"><div class="card-number">{amr_genes}</div><div>AMR Genes</div><i class="fas fa-biohazard fa-2x"></i></div>
        <div class="dashboard-card" onclick="switchTab('virulence')"><div class="card-number">{virulence_genes}</div><div>Virulence Genes</div><i class="fas fa-virus fa-2x"></i></div>
        <div class="dashboard-card" onclick="switchTab('plasmids')"><div class="card-number">{plasmid_genes}</div><div>Plasmid Replicons</div><i class="fas fa-dna fa-2x"></i></div>
        <div class="dashboard-card" onclick="switchTab('patterns')"><div class="card-number">{high_risk_combos}</div><div>High‑Risk Combos</div><i class="fas fa-exclamation-triangle fa-2x"></i></div>
    </div>

    <!-- Tab Navigation -->
    <div class="tab-navigation">
        <button class="tab-button" onclick="switchTab('summary')"><i class="fas fa-chart-pie"></i> Summary</button>
        <button class="tab-button" onclick="switchTab('samples')"><i class="fas fa-list"></i> Samples</button>
        <button class="tab-button" onclick="switchTab('qc')"><i class="fas fa-chart-line"></i> FASTA QC</button>
        <button class="tab-button" onclick="switchTab('mlst')"><i class="fas fa-code-branch"></i> MLST</button>
        <button class="tab-button" onclick="switchTab('amr')"><i class="fas fa-biohazard"></i> AMR</button>
        <button class="tab-button" onclick="switchTab('virulence')"><i class="fas fa-virus"></i> Virulence</button>
        <button class="tab-button" onclick="switchTab('plasmids')"><i class="fas fa-dna"></i> Plasmids</button>
        <button class="tab-button" onclick="switchTab('patterns')"><i class="fas fa-project-diagram"></i> Patterns</button>
        <button class="tab-button" onclick="switchTab('aiguide')"><i class="fas fa-robot"></i> AI Guide</button>
        <button class="tab-button" onclick="switchTab('calltoaction')"><i class="fas fa-globe"></i> Call to Action</button>
        <button class="tab-button" onclick="switchTab('export')"><i class="fas fa-download"></i> Export</button>
    </div>

    <!-- SUMMARY TAB -->
    <div id="summary-tab" class="tab-content active">
        <h2 class="section-header">📊 Executive Summary</h2>
        <div class="info-box">
            <h4><i class="fas fa-info-circle"></i> About this report</h4>
            <p>This gene‑centric report aggregates data from all EnteroMark modules: QC, MLST, AMRfinder, ABRicate (CARD, ResFinder, VFDB, PlasmidFinder, etc.). Each resistance/virulence/plasmid gene is shown with its frequency and the list of genomes that carry it.</p>
            <p><strong>How to use:</strong> Click on any dashboard card or tab to explore detailed data. Use the filters and search bars in each section.</p>
        </div>
        <div class="critical-cards">
            <div class="critical-card" style="border-left-color:#dc2626"><h4><i class="fas fa-skull-crosswalk"></i> Critical Resistance</h4><p>vanA/B/D/M, optrA, cfr, poxtA, high-level aminoglycosides</p></div>
            <div class="critical-card" style="border-left-color:#f59e0b"><h4><i class="fas fa-virus"></i> High‑Risk Virulence</h4><p>esp, ace, gelE, cyl, hyl, acm, scm, biofilm‑associated</p></div>
            <div class="critical-card" style="border-left-color:#10b981"><h4><i class="fas fa-chart-line"></i> QC & Assembly</h4><p>N50, GC%, total bases, etc. – assess assembly quality</p></div>
        </div>
        <div class="alert-box alert-info"><i class="fas fa-lightbulb"></i> <strong>Clinical Impact:</strong> E. faecium is a leading cause of hospital‑acquired infections. Vancomycin‑resistant Enterococcus (VRE) and linezolid‑resistant isolates are critical public health threats. This report helps track resistance genes across isolates.</div>
    </div>

    <!-- SAMPLES TAB (fixed search) -->
    <div id="samples-tab" class="tab-content">
        <h2 class="section-header">📋 Sample Overview</h2>
        <div class="info-box">
            <h4><i class="fas fa-users"></i> Population Structure Overview</h4>
            <p>This table summarises key typing results for each genome. Understanding population structure helps identify dominant clones, track outbreaks, and link genotypes to phenotypes.</p>
            <ul>
                <li><strong>MLST (Sequence Type)</strong>: Gold standard for global epidemiology. Clonal complexes indicate recent common ancestry.</li>
                <li><strong>QC Metrics</strong>: N50, GC%, total bases, etc. – assess assembly completeness and quality.</li>
                <li><strong>ND (Not Determined)</strong>: Indicates that the result could not be determined from available data.</li>
                <li><strong>Click on column headers</strong> to sort the table. Use the search box to filter samples.</li>
            </ul>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="exportTableToCSV('samples-table','samples_overview.csv')"><i class="fas fa-download"></i> Export CSV</button>
            <input type="text" id="search-samples" class="search-box" placeholder="🔍 Search samples..." onkeyup="filterTable('samples-table','search-samples')" style="width: 60%; margin-left: auto;">
        </div>
        <div class="master-scrollable-container">
            <table id="samples-table" class="data-table">
                <thead>
                    <tr>
                        <th>Sample</th><th>ST</th><th>N50</th><th>GC%</th><th>AT%</th><th>N75</th><th>N90</th><th>Total Bases</th>
                    </tr>
                </thead>
                <tbody>{"".join(sample_rows)}</tbody>
            </table>
        </div>
    </div>

    <!-- FASTA QC TAB (fixed sorting) -->
    <div id="qc-tab" class="tab-content">
        <h2 class="section-header">📈 FASTA Quality Control Metrics</h2>
        <div class="info-box">
            <h4><i class="fas fa-chart-line"></i> Assembly Quality Metrics</h4>
            <p>N50, N75, N90: larger values indicate better assembly contiguity. GC% should be around 38-39% for E. faecium. High numbers of short sequences or warnings may indicate poor assembly quality.</p>
            <p><strong>Tip:</strong> Use the search box to find samples with specific QC values (e.g., N50 > 100000). Click any column header to sort.</p>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="exportTableToCSV('qc-table','fasta_qc.csv')"><i class="fas fa-download"></i> Export CSV</button>
            <input type="text" id="search-qc" class="search-box" placeholder="🔍 Search samples..." onkeyup="filterTable('qc-table','search-qc')" style="width: 60%;">
        </div>
        <div class="master-scrollable-container">
            <table id="qc-table" class="data-table">
                <thead>
                    <tr>
                        <th>Sample</th>
                        {''.join([f'<th>{m.replace("_"," ").title()}</th>' for m in qc_headers])}
                    </tr>
                </thead>
                <tbody>{"".join(qc_rows)}</tbody>
            </table>
        </div>
    </div>

    <!-- MLST TAB -->
    <div id="mlst-tab" class="tab-content">
        <h2 class="section-header">🧬 MLST Distribution</h2>
        <div class="info-box">
            <h4><i class="fas fa-code-branch"></i> Sequence Type (ST) Distribution</h4>
            <p>MLST is the gold standard for global epidemiology. Each ST is a unique combination of alleles at seven housekeeping genes. Common STs (e.g., ST17, ST78, ST203) have been associated with hospital outbreaks and different resistance profiles.</p>
            <p><strong>Note:</strong> Click on column headers to sort by frequency or ST number. The "Samples" column lists all genomes belonging to each ST (scrollable).</p>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="exportTableToCSV('mlst-table','mlst.csv')"><i class="fas fa-download"></i> Export CSV</button>
            <input type="text" id="search-mlst" class="search-box" placeholder="🔍 Search ST or sample..." onkeyup="filterTable('mlst-table','search-mlst')" style="width: 60%;">
        </div>
        <div class="master-scrollable-container">
            <table id="mlst-table" class="data-table">
                <thead>
                    <tr><th>ST</th><th>Count</th><th>Percentage</th><th>Samples</th></tr>
                </thead>
                <tbody>{st_rows}</tbody>
            </table>
        </div>
    </div>

    <!-- AMR TAB -->
    <div id="amr-tab" class="tab-content">
        <h2 class="section-header">💊 Antimicrobial Resistance Genes</h2>
        <div class="info-box">
            <h4><i class="fas fa-skull-crosswalk"></i> Clinical Significance of AMR Genes</h4>
            <ul>
                <li><strong>Vancomycin resistance (vanA, vanB, vanD, vanM)</strong>: High‑level resistance to last‑resort glycopeptides. VRE is a WHO high priority pathogen.</li>
                <li><strong>Linezolid resistance (optrA, cfr, poxtA)</strong>: Compromises last‑line oxazolidinones. Often plasmid‑mediated and co‑spreads with vancomycin resistance.</li>
                <li><strong>High‑level aminoglycosides (aac, ant, aph)</strong>: Abolishes synergy with cell‑wall agents, limiting treatment options.</li>
                <li><strong>Beta‑lactamase (pbp5 alterations)</strong>: Intrinsic resistance to ampicillin, common in clinical isolates.</li>
                <li><strong>Other families</strong>: Macrolides (erm), tetracyclines (tet), fluoroquinolones (qnr), etc.</li>
            </ul>
            <p><strong>How to use:</strong> Click category buttons to filter genes. Use the search box to find specific genes, and the genome highlight box to highlight isolates carrying a gene of interest.</p>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="setFilter('van', 'amr-table', 'search-amr')">Vancomycin (van)</button>
            <button class="action-btn" onclick="setFilter('optr', 'amr-table', 'search-amr')">Linezolid (optrA)</button>
            <button class="action-btn" onclick="setFilter('cfr', 'amr-table', 'search-amr')">Linezolid (cfr)</button>
            <button class="action-btn" onclick="setFilter('poxt', 'amr-table', 'search-amr')">Linezolid (poxtA)</button>
            <button class="action-btn" onclick="setFilter('aac', 'amr-table', 'search-amr')">Aminoglycoside</button>
            <button class="action-btn" onclick="setFilter('ant', 'amr-table', 'search-amr')">Aminoglycoside</button>
            <button class="action-btn" onclick="setFilter('aph', 'amr-table', 'search-amr')">Aminoglycoside</button>
            <button class="action-btn" onclick="setFilter('pbp5', 'amr-table', 'search-amr')">Beta‑lactam (pbp5)</button>
            <button class="action-btn" onclick="setFilter('erm', 'amr-table', 'search-amr')">Macrolide</button>
            <button class="action-btn" onclick="setFilter('tet', 'amr-table', 'search-amr')">Tetracycline</button>
            <button class="action-btn" onclick="setFilter('qnr', 'amr-table', 'search-amr')">Fluoroquinolone</button>
            <button class="action-btn" onclick="setFilter('cat', 'amr-table', 'search-amr')">Chloramphenicol</button>
            <button class="action-btn" onclick="setFilter('', 'amr-table', 'search-amr')">Clear</button>
        </div>
        <div class="action-buttons">
            <input type="text" id="search-amr" class="search-box" placeholder="Search genes..." onkeyup="filterTable('amr-table','search-amr')" style="width: 60%;">
            <input type="text" id="highlight-amr" class="search-box" placeholder="Highlight genome (e.g., GCF_001720945.1)" onkeyup="highlightGenomes('amr-table','highlight-amr')" style="width: 40%;">
        </div>
        <div class="master-scrollable-container">
            <table id="amr-table" class="data-table">
                <thead><tr><th>Gene</th><th>Database</th><th>Frequency</th><th>Genomes</th></tr></thead>
                <tbody>{"".join(amr_rows)}</tbody>
            </table>
        </div>
    </div>

    <!-- VIRULENCE TAB -->
    <div id="virulence-tab" class="tab-content">
        <h2 class="section-header">🦠 Virulence Genes</h2>
        <div class="info-box">
            <h4><i class="fas fa-virus"></i> Clinical Significance of Virulence Factors</h4>
            <ul>
                <li><strong>esp (Enterococcal surface protein)</strong>: Promotes biofilm formation and urinary tract colonisation.</li>
                <li><strong>ace (Adhesin of collagen)</strong>: Mediates adherence to extracellular matrix proteins.</li>
                <li><strong>acyl (Cytolysin)</strong>: Pore‑forming toxin that lyses red blood cells and other bacteria.</li>
                <li><strong>gelE (Gelatinase)</strong>: Protease involved in biofilm dispersal and tissue damage.</li>
                <li><strong>acm / scm</strong>: Adhesins to collagen, important for endocarditis pathogenesis.</li>
                <li><strong>fsr locus</strong>: Quorum‑sensing system regulating gelE and sprE.</li>
            </ul>
            <p><strong>Tip:</strong> Use filters to focus on specific virulence mechanisms. The genome highlight feature works across all tabs.</p>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="setFilter('esp', 'virulence-table', 'search-vir')">esp</button>
            <button class="action-btn" onclick="setFilter('ace', 'virulence-table', 'search-vir')">ace</button>
            <button class="action-btn" onclick="setFilter('acm', 'virulence-table', 'search-vir')">acm</button>
            <button class="action-btn" onclick="setFilter('scm', 'virulence-table', 'search-vir')">scm</button>
            <button class="action-btn" onclick="setFilter('cyl', 'virulence-table', 'search-vir')">Cytolysin</button>
            <button class="action-btn" onclick="setFilter('gelE', 'virulence-table', 'search-vir')">gelE</button>
            <button class="action-btn" onclick="setFilter('fsr', 'virulence-table', 'search-vir')">fsr</button>
            <button class="action-btn" onclick="setFilter('hyl', 'virulence-table', 'search-vir')">hyl</button>
            <button class="action-btn" onclick="setFilter('', 'virulence-table', 'search-vir')">Clear</button>
        </div>
        <div class="action-buttons">
            <input type="text" id="search-vir" class="search-box" placeholder="Search virulence genes..." onkeyup="filterTable('virulence-table','search-vir')" style="width: 60%;">
            <input type="text" id="highlight-vir" class="search-box" placeholder="Highlight genome" onkeyup="highlightGenomes('virulence-table','highlight-vir')" style="width: 40%;">
        </div>
        <div class="master-scrollable-container">
            <table id="virulence-table" class="data-table">
                <thead><tr><th>Gene</th><th>Database</th><th>Frequency</th><th>Genomes</th></tr></thead>
                <tbody>{"".join(virulence_rows)}</tbody>
            </table>
        </div>
    </div>

    <!-- PLASMIDS TAB -->
    <div id="plasmids-tab" class="tab-content">
        <h2 class="section-header">🧬 Plasmid Replicons</h2>
        <div class="info-box">
            <h4><i class="fas fa-dna"></i> Plasmid Typing</h4>
            <p>Plasmid replicons indicate the presence of specific plasmid families. Many resistance genes (vanA, optrA, cfr) are carried on mobile plasmids, facilitating horizontal gene transfer.</p>
            <p><strong>Families:</strong> rep (Inc18 family), repA_N (pRE25/pIP501 family), and others. Tracking plasmids helps understand resistance spread.</p>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="setFilter('rep', 'plasmid-table', 'search-plasmid')">Inc18 family</button>
            <button class="action-btn" onclick="setFilter('repA_N', 'plasmid-table', 'search-plasmid')">pRE25 family</button>
            <button class="action-btn" onclick="setFilter('', 'plasmid-table', 'search-plasmid')">Clear</button>
        </div>
        <div class="action-buttons">
            <input type="text" id="search-plasmid" class="search-box" placeholder="Search replicons..." onkeyup="filterTable('plasmid-table','search-plasmid')" style="width: 60%;">
            <input type="text" id="highlight-plasmid" class="search-box" placeholder="Highlight genome" onkeyup="highlightGenomes('plasmid-table','highlight-plasmid')" style="width: 40%;">
        </div>
        <div class="master-scrollable-container">
            <table id="plasmid-table" class="data-table">
                <thead><tr><th>Replicon</th><th>Database</th><th>Frequency</th><th>Genomes</th></tr></thead>
                <tbody>{"".join(plasmid_rows)}</tbody>
            </table>
        </div>
    </div>

    <!-- PATTERNS TAB -->
    <div id="patterns-tab" class="tab-content">
        <h2 class="section-header">🔍 Pattern Discovery</h2>
        <div class="info-box">
            <h4><i class="fas fa-lightbulb"></i> Co‑occurrence Patterns</h4>
            <p>These tables highlight isolates that carry combinations of critical resistance and virulence genes. Such co‑occurrence can indicate high‑risk clones.</p>
            <ul>
                <li><strong>High‑Risk Combinations</strong>: Critical resistance (van, optrA, high‑level aminoglycoside) + high‑risk virulence (esp, gelE, cyl).</li>
                <li><strong>Vancomycin + pbp5 mutations</strong>: May indicate combined glycopeptide and ampicillin resistance.</li>
                <li><strong>Linezolid + Vancomycin resistance</strong>: Very concerning for last‑line treatment failure.</li>
            </ul>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="exportTableToCSV('patterns-table1','high_risk_combos.csv')">Export Table</button>
            <input type="text" id="search-patterns" class="search-box" placeholder="Search patterns..." onkeyup="filterTable('patterns-table1','search-patterns')" style="width: 60%;">
        </div>
        <h3>⚠️ High‑Risk Combinations (Critical Resistance + High‑Risk Virulence)</h3>
        <div class="master-scrollable-container">
            <table id="patterns-table1" class="data-table">
                <thead><tr><th>Sample</th><th>ST</th><th>Critical Resistance</th><th>High‑Risk Virulence</th></tr></thead>
                <tbody>{high_risk_rows}</tbody>
            </table>
        </div>
        <h3>💊 Vancomycin + pbp5 mutations</h3>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Sample</th><th>ST</th><th>van genes</th><th>pbp5 mutations</th></tr></thead><tbody>{van_pbp5_rows}</tbody></table></div>
        <h3>💊 Linezolid resistance + Vancomycin resistance</h3>
        <div class="master-scrollable-container"><table class="data-table"><thead><tr><th>Sample</th><th>ST</th><th>Linezolid genes</th><th>van genes</th></tr></thead><tbody>{linezolid_van_rows}</tbody></table></div>
    </div>

    <!-- AI GUIDE TAB -->
    <div id="aiguide-tab" class="tab-content">
        <h2 class="section-header"><i class="fas fa-robot"></i> AI Assistant Guide</h2>
        <div class="info-box">
            <h4>Using ChatGPT, Claude, or Gemini with this report</h4>
            <p>Upload the <strong>enteromark_ultimate_report.json</strong> file (available in the Export tab) to an AI assistant and ask questions like:</p>
            <ul>
                <li>Which samples carry both vanA and optrA?</li>
                <li>List all isolates with the ST17 plus high‑level aminoglycoside resistance.</li>
                <li>What is the most common virulence gene in this collection?</li>
                <li>Show me the QC metrics for samples with N50 less than 50000.</li>
                <li>Which plasmids are present in VRE isolates?</li>
                <li>Create a summary table of ST distribution and associated resistance genes.</li>
            </ul>
            <p>You can also copy any table from this HTML report into the AI for instant analysis.</p>
        </div>
    </div>

    <!-- CALL TO ACTION TAB -->
    <div id="calltoaction-tab" class="tab-content">
        <h2 class="section-header"><i class="fas fa-globe"></i> Call to Action</h2>
        <div class="info-box" style="background:#fef3c7;">
            <h4>The Global Burden of AMR</h4>
            <p>Antimicrobial resistance (AMR) kills an estimated 1.27 million people annually. <strong>Enterococcus faecium</strong> is a WHO high priority pathogen, with vancomycin‑resistant strains (VRE) causing difficult‑to‑treat hospital infections.</p>
            <p>EnteroMark aims to empower researchers and clinicians to track resistance locally and globally.</p>
        </div>
        <div style="text-align:center; margin:40px 0;">
            <i class="fas fa-star" style="font-size:3em; color:#ffc107;"></i>
            <h3>We invite you to contribute!</h3>
            <p>If you find EnteroMark useful, please <strong>star our GitHub repository</strong> and share it with your network.<br>
            <a href="https://github.com/bbeckley-hub/EnteroMark" target="_blank" style="font-size:1.2em;">https://github.com/bbeckley-hub/EnteroMark</a></p>
            <p>For contributions, bug reports, or feature requests, please open an issue or contact <strong>brownbeckley94@gmail.com</strong>.</p>
            <p><i class="fas fa-chalkboard-user"></i> We welcome collaborations to adapt EnteroMark for other pathogens and to improve AMR surveillance.</p>
        </div>
    </div>

    <!-- EXPORT TAB -->
    <div id="export-tab" class="tab-content">
        <h2 class="section-header">💾 Export Data</h2>
        <div class="info-box">
            <h4>Download data for further analysis</h4>
            <p>You can export each table as CSV or download the complete JSON for AI analysis.</p>
        </div>
        <div class="action-buttons">
            <button class="action-btn" onclick="exportTableToCSV('samples-table','samples.csv')"><i class="fas fa-download"></i> Samples CSV</button>
            <button class="action-btn" onclick="exportTableToCSV('qc-table','fasta_qc.csv')"><i class="fas fa-download"></i> FASTA QC CSV</button>
            <button class="action-btn" onclick="exportTableToCSV('mlst-table','mlst.csv')"><i class="fas fa-download"></i> MLST CSV</button>
            <button class="action-btn" onclick="exportTableToCSV('amr-table','amr.csv')"><i class="fas fa-download"></i> AMR CSV</button>
            <button class="action-btn" onclick="exportTableToCSV('virulence-table','virulence.csv')"><i class="fas fa-download"></i> Virulence CSV</button>
            <button class="action-btn" onclick="exportTableToCSV('plasmid-table','plasmids.csv')"><i class="fas fa-download"></i> Plasmids CSV</button>
            <button class="action-btn" onclick="window.location.href='enteromark_ultimate_report.json'"><i class="fas fa-file-code"></i> Download JSON</button>
        </div>
    </div>

    <!-- COMPREHENSIVE FOOTER -->
    <div class="footer">
        <h3><i class="fas fa-dna"></i> EnteroMark Ultimate Reporter v1.0.0</h3>
        <p>University of Ghana Medical School – Department of Medical Biochemistry</p>
        <p><strong>Author:</strong> Brown Beckley | <strong>Email:</strong> brownbeckley94@gmail.com | <strong>GitHub:</strong> <a href="https://github.com/bbeckley-hub/EnteroMark" target="_blank">https://github.com/bbeckley-hub/EnteroMark</a></p>
        <p><strong>Critical Genes Tracked:</strong> 🔴 Vancomycin (vanA, vanB, vanD, vanM) | 🟠 Linezolid (optrA, cfr, poxtA) | 🟡 High‑level Aminoglycosides (aac, ant, aph) | 🟢 Efflux Pumps & Biocides | 🔵 Adhesins & Biofilm</p>
        <p><i class="fas fa-star"></i> If you find this tool helpful, please star us on GitHub and cite: <em>Beckley et al. (2026) EnteroMark – A comprehensive genomic analysis pipeline for Enterococcus faecium.</em></p>
        <p class="timestamp">Generated on {metadata.get('analysis_date', 'Unknown')}</p>
    </div>
</div>
</body>
</html>"""
        return html


# =============================================================================
# MAIN REPORTER
# =============================================================================
class EnteroMarkUltimateReporter:
    def __init__(self, input_dir: Path):
        self.input_dir = Path(input_dir)
        self.output_dir = self.input_dir / "ENTEROMARK_ULTIMATE_REPORTS"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.parser = EnteroMarkHTMLParser()
        self.analyzer = EnteroMarkDataAnalyzer()
        self.generator = EnteroMarkHTMLGenerator(self.analyzer)
        self.metadata = {
            "tool_name": "EnteroMark Ultimate Reporter",
            "version": "1.0.0",
            "author": "Brown Beckley",
            "affiliation": "University of Ghana Medical School",
            "analysis_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "input_directory": str(self.input_dir)
        }

    def find_files(self) -> Dict[str, Path]:
        files = {}
        qc_candidates = list(self.input_dir.glob("enteromark_fasta_qc_summary.html"))
        if qc_candidates:
            files['qc'] = qc_candidates[0]
        mlst_candidates = list(self.input_dir.glob("mlst_summary.html"))
        if mlst_candidates:
            files['mlst'] = mlst_candidates[0]
        amrfinder_candidates = list(self.input_dir.glob("enteromark_amrfinder_summary_report.html"))
        if amrfinder_candidates:
            files['amrfinder'] = amrfinder_candidates[0]
        files['abricate'] = []
        for db_file in self.parser.db_files:
            candidates = list(self.input_dir.glob(db_file))
            if candidates:
                files['abricate'].append(candidates[0])
        return files

    def run(self):
        print("=" * 80)
        print("🧬 EnteroMark Ultimate Reporter – E. faecium Gene‑Centric Analysis v1.0.0")
        print("=" * 80)
        files = self.find_files()
        if not files:
            print("❌ No required summary HTML files found!")
            return False

        # Parse QC
        qc_data = {}
        if 'qc' in files:
            qc_data = self.parser.parse_qc_summary(files['qc'])
        # Parse MLST
        mlst_data = {}
        if 'mlst' in files:
            mlst_data = self.parser.parse_mlst_summary(files['mlst'])
        # Collect all sample names
        all_samples = set(qc_data.keys()) | set(mlst_data.keys())
        total_samples = len(all_samples)

        # Parse AMRfinder if present
        gene_freqs = {}
        if 'amrfinder' in files:
            amr_freq = self.parser.parse_amrfinder_summary(files['amrfinder'], total_samples)
            if amr_freq:
                gene_freqs['amrfinder'] = amr_freq

        # Parse ABRicate databases
        if 'abricate' in files:
            for db_path in files['abricate']:
                db_freq = self.parser.parse_abricate_database_summary(db_path, total_samples)
                if db_freq:
                    db_name = db_path.stem.replace('enteromark_', '').replace('_summary_report', '')
                    if db_name.endswith('_summary'):
                        db_name = db_name.replace('_summary', '')
                    gene_freqs[db_name] = db_freq

        # Build integrated data
        integrated = {
            'metadata': self.metadata,
            'samples': {s: {} for s in all_samples},
            'qc_data': qc_data,
            'mlst_data': mlst_data,
            'gene_frequencies': gene_freqs,
            'patterns': self.analyzer.create_patterns(
                {s: {'mlst': mlst_data.get(s, {})} for s in all_samples},
                gene_freqs
            )
        }

        # Generate JSON
        json_out = self.output_dir / "enteromark_ultimate_report.json"
        with open(json_out, 'w', encoding='utf-8') as f:
            json.dump(integrated, f, indent=2, default=str)

        # Generate HTML
        self.generator.generate_main_report(integrated, self.output_dir)

        print("✅ All reports generated successfully.")
        return True


def main():
    parser = argparse.ArgumentParser(description='EnteroMark Ultimate Reporter – Gene‑Centric E. faecium Analysis v1.0.0')
    parser.add_argument('-i', '--input-dir', required=True, help='Directory containing EnteroMark summary HTML files')
    args = parser.parse_args()
    reporter = EnteroMarkUltimateReporter(Path(args.input_dir))
    reporter.run()

if __name__ == "__main__":
    main()