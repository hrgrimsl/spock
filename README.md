# SPOCK - Specialized Openfermion Psi4 Chemistry Kit
## To Do
-Take basis set, geometry, frozen orbitals from user

-Get integrals, MO coefficients from Psi4 SCF calculation

-Reorder MO coefficients so that frozen orbitals at top or bottom

-Redo SCF calculation with frozen core/virtuals

-Run FCI calculation with frozen core/virtuals

-Serialize Psi4 data to hdf5

-Tell OFPsi4 to use this Psi4 data for its molecule object

##Feature Requests
-Fix active space issue present in OFPsi4.  We need orbital reordering

-Avoid the process of using psi4 executables.  Writing these instead of using the python API is ugly and inefficient
