import numpy as np

def getMag(v):
  return np.linalg.norm(v)

def getUnitVec(v):
  return v/getMag(v)


