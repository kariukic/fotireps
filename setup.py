#! /usr/bin/env python
"""
Set up for FoTIREPS
"""
from setuptools import setup

requirements = [
    # "scipy>=1.0",
    # others
]

setup(
    name="fotireps",
    version=0.1,
    install_requires=requirements,
    python_requires=">=3.6",
    scripts=[
        "scripts/scripts_writer.py",
        "scripts/rename_di_outputs.py",
    ],
)