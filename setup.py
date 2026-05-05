from setuptools import setup, find_packages

setup(
    name="enteromark",
    version="1.0.0",
    author="Brown Beckley",
    author_email="brownbeckley94@gmail.com",
    description="Complete Enterococcus faecium Genomic Analysis Platform",
    long_description="""EnteroMark – Complete E. faecium genomic analysis pipeline:

🔬 MLST (Multi‑Locus Sequence Typing)
🧬 Antimicrobial Resistance (AMRfinderPlus)
🛡️ Virulence & Resistance Genes (ABRicate with NCBI, CARD, ResFinder, VFDB, etc.)
📊 FASTA Quality Control
📈 Ultimate gene‑centric summary reports

Critical Targets:
🔴 Vancomycin resistance (vanA, vanB, vanD, vanM)
🟠 Linezolid resistance (optrA, cfr, poxtA)
🟡 High‑level aminoglycoside resistance (aac, ant, aph)
🟢 Efflux pumps & biocides (msrC, emeA, qac)
🔵 Adhesins & biofilm formation
""",
    long_description_content_type="text/markdown",
    url="https://github.com/bbeckley-hub/EnteroMark",
    packages=find_packages(include=["enteromark", "enteromark.*"]),
    include_package_data=True,
    package_data={
        "enteromark": [
            "modules/*/*.py",
            "modules/*/bin/*",
            "modules/*/data/**/*",
            "modules/*/db/**/*",
            "modules/*/scripts/*",
            "modules/*/*.fna",         
        ],
    },
    python_requires=">=3.9",
    install_requires=[
        "pandas>=1.5.0",
        "biopython>=1.85",
        "psutil>=5.9.0",
        "requests>=2.28.0",
        "tqdm>=4.64.0",
        "click>=8.0.0",
        "beautifulsoup4>=4.11.0",
        "lxml>=4.9.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.12.0",
    ],
    entry_points={
        "console_scripts": [
            "enteromark=enteromark.enteromark:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ],
)