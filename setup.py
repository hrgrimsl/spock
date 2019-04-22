from distutils.core import setup

setup(
    name='SPOCK',
    version='0.1',
    url='https://github.com/hrgrimsl/spock.git',
    long_description=open('README.txt').read(),
    long_description_content_type='text/markdown',
    author='Harper Grimsley'
    description = 'Specialized Psi4 OpenFermion Chemistry Kit',
    requires_python = '>=3.6.0'
    packages=find_packages(exclude=['*tests*', '*log*']),
)

