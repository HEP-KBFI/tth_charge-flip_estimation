import numpy as np

def calculate_M(Aw):
  M = np.dot(np.linalg.inv( np.dot(np.transpose(Aw), Aw) ) , np.transpose(Aw))
  return M

#Turn ratios list to numpy array
def make_category_matrix(catRatios, weighted = True):
  b = np.array(catRatios)
  if weighted:
    return (b[:,0], b[:,1])
  else:
    return (b[:,0], b[:,1]/b[:,1])

def calculate_rates(M, b):
  return np.dot(M, b)
