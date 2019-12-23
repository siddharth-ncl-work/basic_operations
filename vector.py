import numpy as np

def getMag(v):
  return np.linalg.norm(v)

def getUnitVec(v):
  if getMag(v)==0:
    return [0,0,0]
  return v/getMag(v)

def getDotProduct(v1,v2):
  return np.dot(v1,v2)

def getCrossProduct(v1,v2):
  return np.cross(v1,v2)

def getAngleD(v1,v2):
  theta=math.acos(np.dot(v1,v2)/(getMag(v1)*getMag(v2)))
  return math.degrees(theta)

def getAngleR(v1,v2):
  return np.arccos(np.dot(v1,v2)/(getMag(v1)*getMag(v2)))

def get_dist_pt_line(point,line_vec):
  theta=getAngleR(point,line_vec)
  mag=getMag(point)
  return mag*np.sin(theta)
