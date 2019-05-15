# SPOCK - Specialized Psi4 Openfermion Chemistry Kit
## Installation:  
### Linux/Mac/ARC:
Install with

git clone https://github.com/hrgrimsl/spock.git

Then edit ~/.bashrc to include the line

export PYTHONPATH="${PYTHONPATH}:/home/(path to)/spock"

### Windows:
Use Linux like an adult.

## Usage:
Creates an OpenFermion molecule object with properties according to the following kwargs:

basis- string, basis set to be used; must be a Psi4 standard, e.g. sto-3g, ccpvdz, 6-31g*, etc.  Case-insensitive.

charge- int, charge of the molecule.

geometry- OpenFermion-style geometry nested tuple, e.g. (('H',(0,0,0)),('H',(0,0,1))) in Cartesians, angstroms. 

multiplicity- int, multiplicity of the molecule.  (2S+1 if S is the number of unpaired spins.)  This affects S^2, not necessarily S_z.

active- comma-separated string, e.g. 0,1,2,3,4,5.  Spatial orbitals you want to be active after rotations, etc.  e.g. if you want to swap the 3rd and 6th orbitals, and measure the new 3rd orbital, you include 3 in this string, not 6.

reorder- comma-separated string, original orbitals you want to put in your active space, in the order you want them in.  The mo coefficients corresponding to active will be swapped with these 1 for 1.

output- string, specifies to save the psi4 calculation to {output}.hdf5.  If {output}.hdf5 already exists, that psi4 calculation will be loaded and run instead.

n_fdoccs- Number of frozen occupied orbitals.

occ- comma-separated string, specifies reference state by spin-orbital occupation.  Alphas are even, betas are odd.  This is how you control S_z for a given S^2.

loc- string- True or False, do localize MO's in active space? 

## Data Files Generated
In addition to the actual molecule object, two files will be left in your working directory:

scr.molden- Molden file of the molecular orbitals.  Visualize with jmol.

{output}.hdf5- hdf5 file containing psi4 calculations as well as all necessary input specifications for a molecule.  See OpenFermion for examples of loading from this.


## To Do:
-Automate a test suite

## To Maybe Do:
-Simplify orbital specification

-Add PySCF options
