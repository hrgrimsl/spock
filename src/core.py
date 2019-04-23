import psi4

class molecule:
    def __init__(self, **kwargs):
        print('Initializing molecule.')
        for key in ('geometry', 'basis', 'charge', 'multiplicity', 'active', 'configuration'):
            if key in kwargs:
                setattr(self, key, kwargs[key])
    
    def get_psi_geom(self):
        self.psi_geom = ""
        self.psi_geom += str(self.charge) + ' ' +str(self.multiplicity) + '\n'
        for atom in self.geometry:
            self.psi_geom += atom[0] + ' ' + str(atom[1][0]) + ' ' + str(atom[1][1]) + ' ' + str(atom[1][2]) + '\n'
        print(self.psi_geom)
        
    def scf(self):
        psi4.set_memory('1 GB') 
        psi_molecule = psi4.geometry(self.psi_geom)
        psi4.set_options({'basis': self.basis})
        print(psi4.energy('scf')) 

if __name__ == '__main__':
    GEOMETRY = (('N',(0,0,0)),('N',(0,0,1.5)))
    BASIS = 'STO-3G'
    CHARGE = 0
    MULTIPLICITY = 1
    ACTIVE = None
    CONFIGURATION = None
    molecule = molecule(geometry = GEOMETRY, basis = BASIS, charge = 0, multiplicity = MULTIPLICITY, active = ACTIVE, configuration = CONFIGURATION)
    molecule.get_psi_geom()
    molecule.scf()
