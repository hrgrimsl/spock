from psi4 import core
import os
import shutil
import numpy as np
import copy
import openfermion
from openfermion.config import *
from openfermion.ops import general_basis_change


class molecule:
    def __init__(self, **kwargs):

        for key in ('geometry',
 'basis',
 'charge',
 'multiplicity',
 'active',
 'reorder',
 'n_fdoccs',
 'output',
 'loc',
 'occ'):

            if key in kwargs:
                setattr(self, key, kwargs[key])

        self.get_psi_geom()
        self.molecule = openfermion.hamiltonians.MolecularData(self.geometry, self.basis, self.multiplicity, self.charge)
        self.molecule.active = self.active 
        self.molecule.occ = self.occ

    def get_psi_geom(self):
        self.psi_geom = '\n'
        self.psi_geom += str(self.charge) + ' ' +str(self.multiplicity) + '\n'
        for atom in self.geometry:
            self.psi_geom += atom[0] + ' ' + str(atom[1][0]) + ' ' + str(atom[1][1]) + ' ' + str(atom[1][2]) + '\n'
        self.psi_geom += 'symmetry c1'
    
    def run_psi4(self):
        core.set_output_file(self.output+'.dat')
        import psi4
        psi4.set_memory('1 GB') 
        psi_molecule = psi4.geometry(self.psi_geom)
        psi4.set_options({'basis': self.basis, 'molden_write': False, 'WRITER_FILE_LABEL': str(self.output)+'.dat', 'maxiter': 500, 'fail_on_maxiter': False})
        if self.multiplicity !=1:
            psi4.set_options({'reference': 'ROHF'})

        e, wfn = psi4.energy('scf', return_wfn = True)
        self.molecule.hf_energy = e
        if os.path.exists('./scr.molden'):
            os.system('rm scr.molden')
        if self.active!=None:
            cb = wfn.Cb().to_array()
            ca = wfn.Ca().to_array()
            self.molecule.n_orbitals = len(ca)
            ca[:,self.active+self.reorder]=ca[:,self.reorder+self.active]
            cb[:,self.active+self.reorder]=cb[:,self.reorder+self.active]
            if self.loc == 'True':
                acs = psi4.core.Matrix('null')
                acs = acs.from_array(ca[:,self.active])
                Local = psi4.core.Localizer.build("BOYS", wfn.basisset(), acs)
                Local.localize()
                acs = Local.L
                ca[:,self.active] = acs
                acs2 = psi4.core.Matrix('null2')
                acs2 = acs2.from_array(cb[:,self.active])
                Local = psi4.core.Localizer.build("BOYS", wfn.basisset(), acs2)
                Local.localize()
                acs2 = Local.L
                cb[:,self.active] = acs2
            ca = psi4.core.Matrix.from_array(ca)
            cb = psi4.core.Matrix.from_array(cb)
            wfn.Cb().copy(cb)                     
            wfn.Ca().copy(ca)
            psi4.molden(wfn, 'scr.molden')
            psi4.set_options({'frozen_docc': [self.n_fdoccs]})
            self.n_fnoccs = self.molecule.n_orbitals-self.n_fdoccs-len(self.active)
            psi4.set_options({'frozen_uocc': [self.n_fnoccs]})
            self.CASCI = psi4.energy('fci', ref_wfn = wfn)
            self.molecule.CASCI = self.CASCI
            self.fci_energy = self.CASCI
            self.molecule.fci_energy = self.CASCI
        else:
            self.CASCI = psi4.energy('fci', ref_wfn = wfn)
            self.molecule.CASCI = self.CASCI
            self.fci_energy = self.CASCI
            self.molecule.fci_energy = self.CASCI
            self.molecule.n_orbitals = len(wfn.Ca().to_array())
            psi4.molden(wfn, 'scr.molden')
        self.molecule.nuclear_repulsion = psi_molecule.nuclear_repulsion_energy()
        self.molecule.canonical_orbitals = np.asarray(wfn.Ca())
        if self.active == None:
            self.active = [i for i in range(0, self.molecule.n_orbitals)]
        self.molecule.overlap_integrals = np.asarray(wfn.S())
        self.molecule.n_qubits = 2 * self.molecule.n_orbitals        
        self.molecule.orbital_energies = np.asarray(wfn.epsilon_a())
        self.molecule.fock_matrix = np.asarray(wfn.Fa())
        mints = psi4.core.MintsHelper(wfn.basisset())
        self.molecule.one_body_integrals = general_basis_change(
            np.asarray(mints.ao_kinetic()), self.molecule.canonical_orbitals, (1, 0))
        self.molecule.one_body_integrals += general_basis_change(
            np.asarray(mints.ao_potential()), self.molecule.canonical_orbitals, (1, 0))
        two_body_integrals = np.asarray(mints.ao_eri())
        two_body_integrals.reshape((self.molecule.n_orbitals, self.molecule.n_orbitals,
                                    self.molecule.n_orbitals, self.molecule.n_orbitals))
        two_body_integrals = np.einsum('psqr', two_body_integrals)
        two_body_integrals = general_basis_change(
            two_body_integrals, self.molecule.canonical_orbitals, (1, 1, 0, 0))
        self.molecule.two_body_integrals = two_body_integrals
        doccs = [i for i in range(0, self.n_fdoccs)]

        if self.active!=None:
            self.molecule.n_orbitals = len(self.active)
            self.molecule.n_qubits = 2*self.molecule.n_orbitals
        else:
            self.molecule.n_orbitals = len(self.molecule.canonical_orbitals)
            self.molecule.n_qubits = 2*self.molecule.n_orbitals
        self.molecule.n_fdoccs = self.n_fdoccs
        self.molecule.hamiltonian = self.molecule.get_molecular_hamiltonian(occupied_indices = doccs, active_indices = self.active)
        self.molecule.save()
        return self.molecule

if __name__ == '__main__':
    try:
        shutil.rmtree('scr')
    except:
        pass
    os.mkdir('scr')
    GEOMETRY = [('Sc',(0,0,0))]
    BASIS = 'CC-PVDZ'
    CHARGE = 1
    MULTIPLICITY = 1
    OUTPUT = 'Sc.out'
    ACTIVE = [9,10,11,12,13,14]
    REORDER = [9,10,11,12,13,27]
    N_FDOCCS = 9
    molecule = molecule(geometry = GEOMETRY, basis = BASIS, charge = CHARGE, multiplicity = MULTIPLICITY, active = ACTIVE, reorder = REORDER, n_fdoccs = N_FDOCCS, output = OUTPUT)
    molecule.get_psi_geom()
    molecule.run_psi4()
    

     
