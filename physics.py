import vector
import config
from atomic_mass import atomic_mass_dict

def getMI(atom_cords,atom_type):
  distance=getMag(atom_cords)
  return _getMI([atomic_mass_dict[atom_type]],[distance])  


def _getMI(mass_list,distance_list):
  mr_list=zip(mass_list,distance_list)
  for m,r in mr_list:
    MI+=m*pow(r,2)

  return MI*config.amu*pow(config.angstrom,2)
