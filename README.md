# SPOCK - Specialized Psi4 Openfermion Chemistry Kit
## To Do
-Take basis set, geometry, frozen orbitals from user

-Get integrals, MO coefficients from Psi4 SCF calculation

-Reorder MO coefficients so that frozen orbitals at top or bottom

-Redo SCF calculation with frozen core/virtuals

-Run FCI calculation with frozen core/virtuals

-Serialize Psi4 data to hdf5

-Tell OFPsi4 to use this Psi4 data for its molecule object

-Add PySCF options

##Feature Requests
-Fix active space issue present in OFPsi4.  We need orbital reordering

-Enable orbital localization.  OFPsi4 is awful at this

-Avoid the process of using psi4 executables.  Writing these instead of using the python API is ugly and inefficient
