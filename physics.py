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

def getCom(cords,atom_list=None):
  mass_list=[]
  cords_list=[]
  if atom_list==None:
    atom_list=list(cords['atom_no'].values) 
  for atom_no in atom_list:
    mass=atomic_mass_dict[cords[cords['atom_no']==atom_no]['atom'].values[0]]
    mass_list.append(mass)
    cords_list.append(cords[cords['atom_no']==atom_no][['x','y','z']].values)
  return _getCom(mass_list,cords_list)

def _getCom(mass_list,cords_list):
  com=[0.0,0.0,0.0]
  mr_zip=zip(mass_list,cords_list)
  for m,r in mr_list:
    com[0]+=m*r[0]
    com[1]+=m*r[1]
    com[2]+=m*r[2]
  total_mass=np.sum(mass_list)
  com[0]/=total_mass
  com[1]/=total_mass
  com[2]/=total_mass
  return com

 
