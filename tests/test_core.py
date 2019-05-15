import sys
import os

def test_import():
    sys.path.append('.')
    sys.path.append('..')
    sys.path.append('../..')
    sys.path.append('../../..')
    os.system('ls')
    from spock import core
    molecule = core.molecule(loc = 'False', geometry = (('Li',(0,0,0)),('H',(0,0,1))), basis = 'cc-pvdz', charge = 0, multiplicity = 1, active = [0,1,2,3,4,5], reorder = [0,1,2,3,4,5,6], n_fdoccs = 0, output = 'lihtest', occ = None)  
    molecule.get_psi_geom()
    molecule.run_psi4()
    print('Testing scf...')
    assert(molecule.molecule.hf_energy==-7.881590406211253) 
    print('SCF successful!')
    print('Testing detci...')
    assert(molecule.molecule.fci_energy==-7.881755095606396) 
    print('Detci successful!')
    print('All tests successful!  Consult molden output for a visual sanity check!')
if __name__ == '__main__':
    test_import()
