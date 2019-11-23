import pandas as pd
import numpy as np
import sys
sys.path.extend(['../..'])

from io_chem import io
from atomic_mass import atomic_mass_dict

def getCom(df):
  com={'x':0.0,'y':0.0,'z':0.0}
  df['mass']=df['atom'].apply(str.lower).map(atomic_mass_dict)
  df['com_x']=df['x']*df['mass']
  df['com_y']=df['y']*df['mass']
  df['com_z']=df['z']*df['mass']
  total_mass=df['mass'].sum()
  com['x']=df['com_x'].sum()/total_mass
  com['y']=df['com_y'].sum()/total_mass
  com['z']=df['com_z'].sum()/total_mass
  return com

def _getDist(com,row):
  return np.sqrt((row['x']-com['x'])**2+(row['y']-com['y'])**2+(row['z']-com['z'])**2)   

def getDistFromCom(com,in_df):
  df=in_df.copy()
  dist=[]
  atoms=df.shape[0]
  for i in range(atoms):
    row=df.iloc[i,:][['x','y','z']]
    dist.append(_getDist(com,row))
  df['distance']=dist
  return df

if __name__=='__main__':
  file_path='/home/vanka/siddharth/shailaja_project/al_10_for_center_of_mass'
  df=io.readFile(file_path,file_type='xyz')
  com=getCom(df)
  df=getDistFromCom(com,df)
  print(df)
  print('=================================================================')
  print('COM:\nx = {}\ny = {}\nz = {}'.format(com['x'],com['y'],com['z']))
  print('Atoms at Max distance: {}   Distance:{}'.format(np.argmax(df['distance'].values),df['distance'].max()))
  df.to_csv('output/results_com_dist.csv')


