from setup_tools import setup, find_packages

setup(
    name='SPOCK',
    version='0.1',
    packages = find_packages(),
    long_description=open('README.md').read(),
    author='Harper Grimsley',
    description='Specialized Psi4 OpenFermion Chemistry Kit',
    install_requires = open('requirements.txt', 'r').readlines()
)

