#!/usr/bin/env python3
"""
EnteroMark Orchestrator – Complete E. faecium typing & resistance pipeline
Author: Brown Beckley <brownbeckley94@gmail.com>
Version: 1.0.0 - Fixed ABRicate module, cleaned verbose output
"""

import os
import sys
import glob
import argparse
import subprocess
import shutil
import random
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple, Optional

__version__ = "1.0.0"

# =============================================================================
# Color class (same as before)
# =============================================================================
class Color:
    RESET = '\033[0m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'


class EnteroMarkOrchestrator:
    def __init__(self, verbose: bool = False):
        self.base_dir = Path(__file__).parent.resolve()
        self.modules_dir = self.base_dir / "modules"
        self.verbose = verbose
        self.setup_colors()
        self.quotes = self._get_scientific_quotes()
        self.quote_colors = [
            Color.BRIGHT_CYAN, Color.BRIGHT_GREEN, Color.BRIGHT_YELLOW,
            Color.BRIGHT_MAGENTA, Color.BRIGHT_BLUE, Color.BRIGHT_RED,
            Color.CYAN, Color.GREEN, Color.YELLOW, Color.MAGENTA
        ]
        # Output subdirectory names (within final output folder)
        self.output_dirs = {
            'qc': 'fasta_qc_results',
            'mlst': 'mlst_results',
            'abricate': 'enteromark_abricate_results',
            'amr': 'enteromark_amr_results',
            'summary': 'ENTEROMARK_ULTIMATE_REPORTS'
        }
        # HTML files required for summary (they will be copied to top-level output)
        self.summary_html_files = [
            'enteromark_fasta_qc_summary.html',
            'mlst_summary.html',
            'enteromark_amrfinder_summary_report.html',
            'enteromark_card_summary_report.html',
            'enteromark_ncbi_summary_report.html',
            'enteromark_resfinder_summary_report.html',
            'enteromark_vfdb_summary_report.html',
            'enteromark_argannot_summary_report.html',
            'enteromark_megares_summary_report.html',
            'enteromark_plasmidfinder_summary_report.html',
            'enteromark_bacmet2_summary_report.html'
        ]

    # --------------------------------------------------------------------------
    # Setup colors and printing
    # --------------------------------------------------------------------------
    def setup_colors(self):
        self.color_info = Color.CYAN
        self.color_success = Color.BRIGHT_GREEN
        self.color_warning = Color.BRIGHT_YELLOW
        self.color_error = Color.BRIGHT_RED
        self.color_highlight = Color.BRIGHT_CYAN
        self.color_banner = Color.BRIGHT_MAGENTA
        self.color_module = Color.BRIGHT_BLUE
        self.color_sample = Color.GREEN
        self.color_file = Color.YELLOW
        self.color_reset = Color.RESET

    def print_color(self, message: str, color: str = Color.RESET, bold: bool = False):
        style = Color.BOLD if bold else ''
        print(f"{style}{color}{message}{Color.RESET}")

    def print_header(self, title: str, subtitle: str = ""):
        print()
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_CYAN}{' ' * (40 - len(title)//2)}{title}{Color.RESET}")
        if subtitle:
            print(f"{Color.DIM}{Color.WHITE}{' ' * (42 - len(subtitle)//2)}{subtitle}{Color.RESET}")
        print(f"{Color.BOLD}{Color.BRIGHT_BLUE}{'='*80}{Color.RESET}")
        print()

    def print_info(self, message: str):
        print(f"{self.color_info}[INFO]{Color.RESET} {message}")

    def print_success(self, message: str):
        print(f"{self.color_success}✓{Color.RESET} {message}")

    def print_warning(self, message: str):
        print(f"{self.color_warning}⚠️{Color.RESET} {message}")

    def print_error(self, message: str):
        print(f"{self.color_error}✗{Color.RESET} {message}")

    def print_command(self, cmd_str: str):
        if self.verbose:
            print(f"{Color.DIM}{Color.WHITE}  $ {cmd_str}{Color.RESET}")

    # --------------------------------------------------------------------------
    # Quotes (unchanged)
    # --------------------------------------------------------------------------
    def _get_scientific_quotes(self):
        return [
            {"quote": "Science is organised knowledge.", "author": "Herbert Spencer"},
            {"quote": "Biology is the most powerful technology ever created.", "author": "Freeman Dyson"},
            {"quote": "Genomics is a lens on biology.", "author": "Eric Lander"},
            {"quote": "Every microbe has its own story.", "author": "Anonymous"},
            {"quote": "Data beats emotions.", "author": "Sean Rad"},
            {"quote": "Microbes rule the world.", "author": "Paul Stamets"},
            {"quote": "Genes are the language of life.", "author": "Francis Collins"},
            {"quote": "Evolution in a petri dish.", "author": "Richard Lenski"},
            {"quote": "The good thing about science is that it's true whether or not you believe in it.", "author": "Neil deGrasse Tyson"},
            {"quote": "Nothing in life is to be feared, it is only to be understood.", "author": "Marie Curie"},
        ]

    def display_random_quote(self):
        if not self.quotes:
            return
        quote_data = random.choice(self.quotes)
        quote = quote_data["quote"]
        author = quote_data["author"]
        quote_color = random.choice(self.quote_colors)
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print()
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}[{current_time}] 💡 SCIENTIFIC INSIGHT:{Color.RESET}")
        print()
        print(f"{quote_color}   \"{quote}\"{Color.RESET}")
        print(f"{Color.BOLD}{Color.WHITE}   — {author}{Color.RESET}")
        print(f"{Color.DIM}{Color.WHITE}{'─' * 80}{Color.RESET}")
        print()

    # --------------------------------------------------------------------------
    # File discovery helpers
    # --------------------------------------------------------------------------
    def find_fasta_files(self, input_path: str) -> List[Path]:
        self.print_info(f"Searching for files with pattern: {input_path}")
        if '*' in input_path or '?' in input_path:
            matched_files = glob.glob(input_path)
            fasta_files = [Path(f) for f in matched_files if Path(f).is_file() and
                           f.lower().endswith(('.fna', '.fasta', '.fa', '.fn')) and
                           not Path(f).name.startswith('.')]
            self.print_success(f"Found {len(fasta_files)} FASTA files")
            return sorted(fasta_files)
        input_path_obj = Path(input_path)
        if input_path_obj.is_file() and input_path_obj.suffix.lower() in ['.fna', '.fasta', '.fa', '.fn']:
            self.print_success(f"Found single FASTA file: {input_path_obj.name}")
            return [input_path_obj]
        if input_path_obj.is_dir():
            patterns = [f"{input_path}/*.fna", f"{input_path}/*.fasta", f"{input_path}/*.fa", f"{input_path}/*.fn"]
            fasta_files = []
            for pattern in patterns:
                matched_files = glob.glob(pattern)
                for file_path in matched_files:
                    path = Path(file_path)
                    if path.is_file() and not path.name.startswith('.'):
                        fasta_files.append(path)
            fasta_files = sorted(list(set(fasta_files)))
            if fasta_files:
                self.print_success(f"Found {len(fasta_files)} FASTA files in directory")
            else:
                self.print_warning(f"No FASTA files found in directory: {input_path}")
            return fasta_files
        self.print_error(f"Input path not found: {input_path}")
        return []

    def get_file_pattern(self, fasta_files: List[Path]) -> str:
        if not fasta_files:
            return "*.fna"
        extensions = set(f.suffix.lower() for f in fasta_files)
        if len(extensions) == 1:
            return f"*{list(extensions)[0]}"
        return "*"

    # --------------------------------------------------------------------------
    # AMR Database Management (quiet version)
    # --------------------------------------------------------------------------
    def update_amr_database(self) -> bool:
        amr_module_path = self.modules_dir / "amr_module"
        amr_script = amr_module_path / "enteromark_amr.py"
        if not amr_script.exists():
            self.print_error(f"AMR script not found at: {amr_script}")
            return False
        self.print_info("Updating AMRfinderPlus database...")
        cmd = [sys.executable, str(amr_script), "--update-db"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0:
            self.print_success("AMR database updated successfully.")
            return True
        else:
            self.print_error("AMR database update failed.")
            if result.stderr:
                self.print_error(result.stderr.strip())
            return False

    def ensure_amr_database(self) -> bool:
        amr_module_path = self.modules_dir / "amr_module"
        amr_script = amr_module_path / "enteromark_amr.py"
        if not amr_script.exists():
            self.print_error("AMR script not found, cannot check database.")
            return False
        # Query version quietly – we'll parse only the version number
        cmd = [sys.executable, str(amr_script), "--db-version"]
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
        if result.returncode == 0 and "Unknown" not in result.stdout and "No database" not in result.stdout:
            version = result.stdout.strip().split('\n')[-1]  # take last line
            self.print_success(f"AMR database ready (version: {version})")
            return True
        else:
            self.print_warning("AMR database not found or outdated. Attempting automatic update...")
            return self.update_amr_database()

    # --------------------------------------------------------------------------
    # Cleanup module directory (remove copied FASTA files and temporary outputs)
    # --------------------------------------------------------------------------
    def cleanup_module(self, module_path: Path, fasta_files: List[Path]):
        try:
            for fasta_file in fasta_files:
                temp = module_path / fasta_file.name
                if temp.exists():
                    temp.unlink()
            # Remove output subdirectories that may have been created
            for subdir in self.output_dirs.values():
                dir_path = module_path / subdir
                if dir_path.exists():
                    shutil.rmtree(dir_path)
            # Also remove any stray HTML files
            for html in module_path.glob("*.html"):
                html.unlink()
        except Exception as e:
            if self.verbose:
                self.print_warning(f"Cleanup issue in {module_path.name}: {e}")

    # --------------------------------------------------------------------------
    # Module runners (dedicated, following AcinetoScope pattern)
    # --------------------------------------------------------------------------
    def run_qc(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        qc_module_path = self.modules_dir / "fastaqc_module"
        log = ""
        try:
            self.print_header("FASTA QC ANALYSIS", "Quality Control")
            qc_script = qc_module_path / "enteromark_fasta_qc.py"
            if not qc_script.exists():
                return False, f"QC script missing: {qc_script}"
            # Copy FASTA files
            for f in fasta_files:
                shutil.copy2(f, qc_module_path / f.name)
            log += f"Copied {len(fasta_files)} files to QC module\n"
            pattern = self.get_file_pattern(fasta_files)
            pattern_clean = pattern.strip('"')
            cmd = [sys.executable, str(qc_script), pattern_clean]
            self.print_command(' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=qc_module_path)
            if result.returncode != 0:
                log += f"⚠️ QC had warnings (exit {result.returncode})\n"
                if result.stderr:
                    log += result.stderr + "\n"
            # Copy results
            source = qc_module_path / self.output_dirs['qc']
            target = output_dir / self.output_dirs['qc']
            if source.exists():
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(source, target)
                log += f"✓ QC results → {target}\n"
                # Copy main HTML report to top-level output
                html_src = source / "enteromark_fasta_qc_summary.html"
                if html_src.exists():
                    shutil.copy2(html_src, output_dir / "enteromark_fasta_qc_summary.html")
                    log += "✓ Copied QC HTML report\n"
                else:
                    log += "⚠️ QC HTML report missing\n"
            else:
                log += f"⚠️ QC results not found at {source}\n"
            self.display_random_quote()
            return True, log
        except Exception as e:
            return False, f"QC exception: {str(e)}"
        finally:
            self.cleanup_module(qc_module_path, fasta_files)

    def run_mlst(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        mlst_module_path = self.modules_dir / "mlst_module"
        log = ""
        try:
            self.print_header("MLST ANALYSIS", "Multi-Locus Sequence Typing")
            mlst_script = mlst_module_path / "enteromark_mlst.py"
            if not mlst_script.exists():
                return False, f"MLST script missing: {mlst_script}"
            for f in fasta_files:
                shutil.copy2(f, mlst_module_path / f.name)
            log += f"Copied {len(fasta_files)} files to MLST module\n"
            pattern = self.get_file_pattern(fasta_files)
            pattern_clean = pattern.strip('"')
            db_path = mlst_module_path / "db" / "pubmlst"
            out_subdir = self.output_dirs['mlst']
            cmd = [sys.executable, str(mlst_script), "-i", pattern_clean, "-o", out_subdir,
                   "-db", str(db_path), "-sc", str(mlst_module_path), "--batch"]
            self.print_command(' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=mlst_module_path)
            if result.returncode != 0:
                log += f"⚠️ MLST had warnings (exit {result.returncode})\n"
                if result.stderr:
                    log += result.stderr + "\n"
            source = mlst_module_path / out_subdir
            target = output_dir / out_subdir
            if source.exists():
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(source, target)
                log += f"✓ MLST results → {target}\n"
                html_src = source / "mlst_summary.html"
                if html_src.exists():
                    shutil.copy2(html_src, output_dir / "mlst_summary.html")
                    log += "✓ Copied MLST HTML report\n"
                else:
                    log += "⚠️ MLST HTML report missing\n"
            else:
                log += f"⚠️ MLST results not found at {source}\n"
            self.display_random_quote()
            return True, log
        except Exception as e:
            return False, f"MLST exception: {str(e)}"
        finally:
            self.cleanup_module(mlst_module_path, fasta_files)

    def run_abricate(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        ab_module_path = self.modules_dir / "abricate_module"
        log = ""
        try:
            self.print_header("ABRICATE ANALYSIS", "Resistance & Virulence Screening")
            ab_script = ab_module_path / "enteromark_abricate.py"
            if not ab_script.exists():
                return False, f"ABRicate script missing: {ab_script}"
            for f in fasta_files:
                shutil.copy2(f, ab_module_path / f.name)
            log += f"Copied {len(fasta_files)} files to ABRicate module\n"
            pattern = self.get_file_pattern(fasta_files)
            pattern_clean = pattern.strip('"')
            out_subdir = self.output_dirs['abricate']
            # Run the script with explicit output directory
            cmd = [sys.executable, str(ab_script), pattern_clean, "-o", out_subdir]
            self.print_command(' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=ab_module_path)
            if result.returncode != 0:
                log += f"⚠️ ABRicate had warnings (exit {result.returncode})\n"
                if result.stderr:
                    log += result.stderr + "\n"
            source = ab_module_path / out_subdir
            target = output_dir / out_subdir
            if source.exists():
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(source, target)
                log += f"✓ ABRicate results → {target}\n"
                # Copy each expected HTML report to top-level output
                copied = 0
                for html in self.summary_html_files:
                    if html.startswith('enteromark_') and 'card' in html or 'ncbi' in html or 'resfinder' in html or 'vfdb' in html or 'argannot' in html or 'megares' in html or 'plasmidfinder' in html or 'bacmet2' in html:
                        src_file = source / html
                        if src_file.exists():
                            shutil.copy2(src_file, output_dir / html)
                            copied += 1
                log += f"Copied {copied} ABRicate HTML reports\n"
            else:
                log += f"⚠️ ABRicate results not found at {source}\n"
            self.display_random_quote()
            return True, log
        except Exception as e:
            return False, f"ABRicate exception: {str(e)}"
        finally:
            self.cleanup_module(ab_module_path, fasta_files)

    def run_amr(self, fasta_files: List[Path], output_dir: Path, threads: int) -> Tuple[bool, str]:
        if not self.ensure_amr_database():
            return False, "AMR database missing and update failed."
        amr_module_path = self.modules_dir / "amr_module"
        log = ""
        try:
            self.print_header("AMR ANALYSIS", "AMRFinderPlus Resistance Detection")
            amr_script = amr_module_path / "enteromark_amr.py"
            if not amr_script.exists():
                return False, f"AMR script missing: {amr_script}"
            for f in fasta_files:
                shutil.copy2(f, amr_module_path / f.name)
            log += f"Copied {len(fasta_files)} files to AMR module\n"
            pattern = self.get_file_pattern(fasta_files)
            pattern_clean = pattern.strip('"')
            out_subdir = self.output_dirs['amr']
            cmd = [sys.executable, str(amr_script), pattern_clean, "-o", out_subdir]
            self.print_command(' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=amr_module_path)
            if result.returncode != 0:
                log += f"⚠️ AMR had warnings (exit {result.returncode})\n"
                if result.stderr:
                    log += result.stderr + "\n"
            source = amr_module_path / out_subdir
            target = output_dir / out_subdir
            if source.exists():
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(source, target)
                log += f"✓ AMR results → {target}\n"
                html_src = source / "enteromark_amrfinder_summary_report.html"
                if html_src.exists():
                    shutil.copy2(html_src, output_dir / "enteromark_amrfinder_summary_report.html")
                    log += "✓ Copied AMR HTML report\n"
                else:
                    log += "⚠️ AMR HTML report missing\n"
            else:
                log += f"⚠️ AMR results not found at {source}\n"
            self.display_random_quote()
            return True, log
        except Exception as e:
            return False, f"AMR exception: {str(e)}"
        finally:
            self.cleanup_module(amr_module_path, fasta_files)

    def run_summary(self, output_dir: Path) -> Tuple[bool, str]:
        summary_module_path = self.modules_dir / "summary_module"
        log = ""
        try:
            self.print_header("ULTIMATE REPORTER", "Gene‑centric Integrated Analysis")
            summary_script = summary_module_path / "enteromark_summary.py"
            if not summary_script.exists():
                return False, f"Summary script missing: {summary_script}"
            # Copy required HTML files (all must be present in output_dir now)
            copied = 0
            missing = []
            for html in self.summary_html_files:
                src = output_dir / html
                if src.exists():
                    shutil.copy2(src, summary_module_path / html)
                    copied += 1
                else:
                    missing.append(html)
            if missing:
                log += f"⚠️ Missing {len(missing)} HTML files: {', '.join(missing)}\n"
            else:
                log += f"✓ Copied {copied} HTML files to summary module\n"
            # Run the summary reporter
            cmd = [sys.executable, str(summary_script), "-i", "."]
            self.print_command(' '.join(cmd))
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=summary_module_path)
            log += result.stdout
            if result.stderr:
                log += result.stderr
            if result.returncode != 0:
                log += f"⚠️ Summary reporter had issues (exit {result.returncode})\n"
            # Copy ultimate reports
            source = summary_module_path / self.output_dirs['summary']
            target = output_dir / self.output_dirs['summary']
            if source.exists():
                if target.exists():
                    shutil.rmtree(target)
                shutil.copytree(source, target)
                log += f"✓ Ultimate reports → {target}\n"
            else:
                log += f"⚠️ Ultimate reports not found at {source}\n"
            self.display_random_quote()
            return result.returncode == 0, log
        except Exception as e:
            return False, f"Summary exception: {str(e)}"

    # --------------------------------------------------------------------------
    # Clean all modules after run
    # --------------------------------------------------------------------------
    def cleanup_all_modules(self, fasta_files: List[Path]):
        for module_dir in ["fastaqc_module", "mlst_module", "abricate_module", "amr_module"]:
            module_path = self.modules_dir / module_dir
            if module_path.exists():
                self.cleanup_module(module_path, fasta_files)

    # --------------------------------------------------------------------------
    # Main orchestration
    # --------------------------------------------------------------------------
    def run_complete_analysis(self, input_path: str, output_dir: str, threads: int = 1,
                              skip_modules: Dict[str, bool] = None, skip_summary: bool = False,
                              update_amr_db_only: bool = False):
        if update_amr_db_only:
            self.update_amr_database()
            return
        if skip_modules is None:
            skip_modules = {}
        start_time = datetime.now()
        self.display_banner()
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        fasta_files = self.find_fasta_files(input_path)
        if not fasta_files:
            self.print_error("No FASTA files found! Analysis stopped.")
            return
        self.print_success(f"Starting analysis of {len(fasta_files)} E. faecium samples")
        # Create output subdirectories
        for subdir in self.output_dirs.values():
            (output_path / subdir).mkdir(exist_ok=True)
        # Display plan
        self.print_header("ANALYSIS PLAN", "Modules to be executed")
        plan = [
            ("FASTA QC", not skip_modules.get('qc', False)),
            ("MLST", not skip_modules.get('mlst', False)),
            ("ABRicate", not skip_modules.get('abricate', False)),
            ("AMRfinder", not skip_modules.get('amr', False)),
            ("Ultimate Reporter", not skip_summary),
        ]
        for analysis, enabled in plan:
            if enabled:
                print(f"   {Color.BRIGHT_GREEN}✅ ENABLED{Color.RESET} - {analysis}")
            else:
                print(f"   {Color.YELLOW}⏸️  SKIPPED{Color.RESET} - {analysis}")
        print()
        # Run analyses sequentially (to avoid conflicts and keep logs clean)
        results = {}
        if not skip_modules.get('qc', False):
            success, log = self.run_qc(fasta_files, output_path, threads)
            results['QC'] = success
            print(log)
        if not skip_modules.get('mlst', False):
            success, log = self.run_mlst(fasta_files, output_path, threads)
            results['MLST'] = success
            print(log)
        if not skip_modules.get('abricate', False):
            success, log = self.run_abricate(fasta_files, output_path, threads)
            results['ABRicate'] = success
            print(log)
        if not skip_modules.get('amr', False):
            success, log = self.run_amr(fasta_files, output_path, threads)
            results['AMR'] = success
            print(log)
        if not skip_summary:
            success, log = self.run_summary(output_path)
            results['Summary'] = success
            print(log)
        # Clean up module directories
        self.cleanup_all_modules(fasta_files)
        # Final summary
        analysis_time = datetime.now() - start_time
        self.print_header("ANALYSIS COMPLETE", f"Time elapsed: {str(analysis_time).split('.')[0]}")
        self.print_success(f"🎉 Analysis complete! Results in: {output_path}")
        for subdir in sorted(output_path.iterdir()):
            if subdir.is_dir():
                file_count = len(list(subdir.glob("*")))
                self.print_info(f"  📁 {subdir.name} ({file_count} files)")
        self.display_random_quote()

    # --------------------------------------------------------------------------
    # Banner and help
    # --------------------------------------------------------------------------
    def display_banner(self):
        banner = f"""{Color.BOLD}{Color.BRIGHT_MAGENTA}
{'='*80}
{' '*22}🦠 ENTEROMARK - E. faecium Genomic Analysis Pipeline v{__version__}
{'='*80}
{Color.RESET}{Color.BRIGHT_CYAN}
Complete Enterococcus faecium genomic analysis pipeline
MLST | AMR | Virulence | Plasmids | Quality Control | Summary Reports

Critical Genes Tracked:
🔴 Vancomycin (vanA, vanB, vanD, vanM)
🟠 Linezolid (optrA, cfr, poxtA)
🟡 High-level Aminoglycosides (aac, ant, aph)
🟢 Efflux Pumps & Biocides
🔵 Adhesins & Biofilm
{Color.RESET}{Color.DIM}
Author: Brown Beckley | Email: brownbeckley94@gmail.com
Affiliation: University of Ghana Medical School - Department of Medical Biochemistry
Version: {__version__}
{'='*80}{Color.RESET}
"""
        print(banner)

    def print_colored_help(self):
        self.display_banner()
        print(f"{Color.BRIGHT_YELLOW}USAGE:{Color.RESET}")
        print(f"  {Color.GREEN}enteromark{Color.RESET} {Color.CYAN}-i INPUT -o OUTPUT{Color.RESET} [OPTIONS]")
        print(f"  {Color.GREEN}enteromark --update-amr-db{Color.RESET}")
        print()
        print(f"{Color.BRIGHT_YELLOW}REQUIRED ARGUMENTS:{Color.RESET}")
        print(f"  {Color.GREEN}-i, --input{Color.RESET} INPUT    FASTA file(s) (glob pattern like \"*.fna\")")
        print(f"  {Color.GREEN}-o, --output{Color.RESET} OUTPUT  Output directory\n")
        print(f"{Color.BRIGHT_YELLOW}OPTIONAL:{Color.RESET}")
        print(f"  {Color.GREEN}-t, --threads{Color.RESET} THREADS  Threads (default: 2)")
        print(f"  {Color.GREEN}--skip-qc{Color.RESET}             Skip QC")
        print(f"  {Color.GREEN}--skip-mlst{Color.RESET}           Skip MLST")
        print(f"  {Color.GREEN}--skip-abricate{Color.RESET}       Skip ABRicate")
        print(f"  {Color.GREEN}--skip-amr{Color.RESET}            Skip AMR")
        print(f"  {Color.GREEN}--skip-summary{Color.RESET}        Skip final report")
        print(f"  {Color.GREEN}--update-amr-db{Color.RESET}       Update AMR database only")
        print(f"  {Color.GREEN}--verbose{Color.RESET}             Show full commands")
        print(f"  {Color.GREEN}--version{Color.RESET}             Show version")
        print(f"  {Color.GREEN}-h, --help{Color.RESET}            This help\n")
        print(f"{Color.BRIGHT_YELLOW}EXAMPLES:{Color.RESET}")
        print(f"  {Color.GREEN}enteromark -i \"*.fna\" -o results{Color.RESET}")
        print(f"  {Color.GREEN}enteromark -i genome.fna -o results --skip-abricate{Color.RESET}")


# =============================================================================
# Main entry point
# =============================================================================
def main():
    if '--version' in sys.argv:
        print(f"enteromark {__version__}")
        sys.exit(0)
    if '-h' in sys.argv or '--help' in sys.argv:
        orchestrator = EnteroMarkOrchestrator(verbose=False)
        orchestrator.print_colored_help()
        sys.exit(0)
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-t', '--threads', type=int, default=2)
    parser.add_argument('--skip-qc', action='store_true')
    parser.add_argument('--skip-mlst', action='store_true')
    parser.add_argument('--skip-abricate', action='store_true')
    parser.add_argument('--skip-amr', action='store_true')
    parser.add_argument('--skip-summary', action='store_true')
    parser.add_argument('--update-amr-db', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    args = parser.parse_args()
    if args.update_amr_db:
        orchestrator = EnteroMarkOrchestrator(verbose=args.verbose)
        orchestrator.run_complete_analysis("", "", update_amr_db_only=True)
        sys.exit(0)
    if not args.input or not args.output:
        print(f"{Color.BRIGHT_RED}Error: -i and -o are required for analysis.{Color.RESET}")
        sys.exit(1)
    skip_modules = {
        'qc': args.skip_qc,
        'mlst': args.skip_mlst,
        'abricate': args.skip_abricate,
        'amr': args.skip_amr,
    }
    orchestrator = EnteroMarkOrchestrator(verbose=args.verbose)
    orchestrator.run_complete_analysis(
        input_path=args.input,
        output_dir=args.output,
        threads=args.threads,
        skip_modules=skip_modules,
        skip_summary=args.skip_summary
    )

if __name__ == "__main__":
    main()