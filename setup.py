from setuptools import setup, find_packages
import os
import sys


setup(
    name='SPOCK',
    version='0.1',
    packages = find_packages(),
    long_description=open('README.md').read(),
    author='Harper Grimsley',
    author_email='hrgrimsl@vt.edu',
    url='https://github.com/hrgrimsl/spock.git',
    description='Specialized Psi4 OpenFermion Chemistry Kit',
    install_requires = open('requirements.txt', 'r').readlines()
)



