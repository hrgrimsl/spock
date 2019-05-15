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

active- list, e.g. [0,1,2,3,4,5].  Spatial orbitals you want to be active after rotations, etc.  e.g. if you want to swap the 3rd and 6th orbitals, and measure the new 3rd orbital, you include 3 in this string, not 6.

reorder- list, original orbitals you want to put in your active space, in the order you want them in.  The mo coefficients corresponding to active will be swapped with these 1 for 1.

output- string, specifies the .dat psi4 file to be written to, as well as the molden file.  

n_fdoccs- Number of frozen occupied orbitals.

occ- list, specifies reference state by spin-orbital occupation.  Alphas are even, betas are odd.  This is how you control S_z for a given S^2.

loc- string- True or False, do localize MO's in active space? 

## Data Files Generated
In addition to the actual molecule object, two files will be left in your working directory:

{output}.molden - Molden file of the molecular orbitals.  Visualize with jmol.

{output}.dat- Psi4 output file.

## Thoughts:
-Simplify orbital specification

-Add PySCF options

## See Also:

OpenFermionPsi4 - A very similar code to this which I felt was not generally appropriate for chemists, hence the design of this code.  Given their time investment and resources, their program is much nicer than this if you're not trying to do something outside of its wheelhouse.

OpenFermion - I don't know why you're using this if you don't know what OpenFermion is.

Psi4 - Quantum chemistry program which this acts as a wrapper for.  Free, god-tier API, and versatile.  You should donate to its devs!

JMOL - Quantum chemistry visualizer in the vein of Avagadro, VMD, etc.  Useful for visual checks that your mo's are what you want them to be.


