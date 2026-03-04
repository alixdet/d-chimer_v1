#!/usr/bin/env python3
"""
Setup script for d-chimer package.

d-chimer is a pipeline to detect and assign fragments of contigs included
in chimeric sequences during taxonomic assignment with BLAST.
"""

from setuptools import setup, find_packages
import os

# Read the README file
with open('README.md', 'r', encoding='utf-8') as f:
    long_description = f.read()

# Read requirements from environment.yml-equivalent
install_requires = [
    'biopython>=1.81',
    'PyYAML>=5.0',
]

setup(
    name='dchimer',
    version='1.0.0',
    description='Detect and assign chimeric sequences in taxonomic assignment with BLAST',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Sourakhata Tirera, Jean-Marc Frigerio, Alix de Thoisy',
    author_email='alixdet@protonmail.com',
    license='GNU General Public License v3.0',
    url='https://github.com/yourusername/d-chimer',  # Update with actual repo URL
    packages=find_packages(exclude=['tests', '__pycache__']),
    include_package_data=True,
    package_data={
        'dchimer': [
            'dchimer_config.yaml',
            'add_taxo.sh',
            'taxonomy/*.dmp',
        ],
    },
    python_requires='>=3.7',
    install_requires=install_requires,
    entry_points={
        'console_scripts': [
            'dchimer=dchimer.dchimer:main',
        ],
    },
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    keywords='bioinformatics, blast, taxonomy, chimera detection',
    project_urls={
        'Documentation': 'https://github.com/yourusername/d-chimer#readme',
        'Source': 'https://github.com/yourusername/d-chimer',
        'Tracker': 'https://github.com/yourusername/d-chimer/issues',
    },
)
