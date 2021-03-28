#!/usr/bin/env python

from setuptools import setup, find_packages
with open("README.md", 'r') as f:
    long_description = f.read()
    
setup(
    name='titrationFitter',
    version='0.1.2',
    python_requires='>3.7.0',
    description='Extracts binding constant from titration experiments',
    author='Frédéric Dux',
    author_email="duxfrederic@gmail.com",
    url="https://gitlab.epfl.ch/dux/titration-fitting/",
    long_description=long_description,
    license='MIT',
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib"
    ]
)
