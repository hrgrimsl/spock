# SPOCK - Specialized Psi4 Openfermion Chemistry Kit
## Installation:  
### Linux/Mac/ARC:
Install with

git clone https://github.com/hrgrimsl/spock.git

Then edit ~/.bashrc to include the line

export PYTHONPATH="${PYTHONPATH}:/home/hrgrimsl/spock"

### Windows:
Use Linux like an adult.

## Usage:
Creates an OpenFermion molecule object with properties according to the following kwargs:

basis- string, basis set to be used; must be a Psi4 standard, e.g. sto-3g, ccpvdz, 6-31g*, etc.  Case-insensitive.

charge- int, charge of the molecule.

multiplicity- int, multiplicity of the molecule.  (2S+1 if S is the number of unpaired spins.)  This affects S^2, not necessarily S_z.

active- comma-separated string, e.g. 0,1,2,3,4,5.  Spatial orbitals you want to be active after rotations, etc.  e.g. if you want to swap the 3rd and 6th orbitals, and measure the new 3rd orbital, you include 3 in this string, not 6.

reorder- comma-separated string, original orbitals you want to put in your active space, in the order you want them in.  The mo coefficients corresponding to active will be swapped with these 1 for 1.

occ- comma-separated string, specifies reference state by spin-orbital occupation.  Alphas are even, betas are odd.  This is how you control S_z for a given S^2.

loc- string- True or False, do localize MO's in active space? 


## To Do:
-Add PySCF options

-Automate a test suite

-Automate installation

-Introduce better control over saving more than the most recent intermediate calculations.

-Automate cleanup of old files.

-Spice up the log file to be more verbose and not just a data column?

-Simplify orbital specification?

-Remove deprecated kwargs
