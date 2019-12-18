import numpy as np

from . import vector
from . import atomic_mass
from . import constants

def getMI(atom_cords,atom_type,axis):
  distance=vector.getMag(atom_cords)
  return _getMI([atom_cords],[atomic_mass_dict[atom_type.lower()]],axis)  

def _getMI(cords_list,mass_list,axis):
  MI=0
  mr_list=zip(mass_list,cords_list)
  for m,cords in mr_list:
    r=vector.get_dist_pt_line(cords,axis)
    MI+=m*pow(r,2)
  return MI*config.amu*pow(config.angstrom,2)

def getCom(cords,atom_list=None):
  mass_list=[]
  cords_list=[]
  if atom_list==None:
    atom_list=list(cords['atom_no'].values) 
  for atom_no in atom_list:
    mass=atomic_mass.atomic_mass_dict[cords[cords['atom_no']==atom_no]['atom'].values[0].lower()]
    mass_list.append(mass)
    cords_list.append(cords[cords['atom_no']==atom_no][['x','y','z']].values[0])
  return _getCom(cords_list,mass_list)

def _getCom(cords_list,mass_list):
  com=[0.0,0.0,0.0]
  mr_list=zip(mass_list,cords_list)
  for m,r in mr_list:
    com[0]+=m*r[0]
    com[1]+=m*r[1]
    com[2]+=m*r[2]
  total_mass=np.sum(mass_list)
  com[0]/=total_mass
  com[1]/=total_mass
  com[2]/=total_mass
  return com

def getTotalMass(cords,atom_list=None):
  total_mass=0
  if atom_list==None:
    atom_list=list(cords['atom_no'].values)
  for atom_no in atom_list:
    mass=atomic_mass.atomic_mass_dict[cords[cords['atom_no']==atom_no]['atom'].values[0].lower()]
    total_mass+=mass
  return total_mass*constants.amu

