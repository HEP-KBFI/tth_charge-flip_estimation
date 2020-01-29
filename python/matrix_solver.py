import numpy as np
import math

from tthAnalysis.ChargeFlipEstimation.utils import fit_results_to_file

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

def make_coefficient_matrix(exclude_bins = None):
 coeffs = [
   [ 2,  0,  0,  0,  0,  0],
   [ 1,  1,  0,  0,  0,  0],
   [ 0,  2,  0,  0,  0,  0],
   [ 1,  0,  1,  0,  0,  0],
   [ 0,  1,  1,  0,  0,  0],
   [ 0,  0,  2,  0,  0,  0],
   [ 0,  0,  0,  2,  0,  0],
   [ 0,  0,  0,  1,  1,  0],
   [ 0,  0,  0,  0,  2,  0],
   [ 0,  0,  0,  1,  0,  1],
   [ 0,  0,  0,  0,  1,  1],
   [ 0,  0,  0,  0,  0,  2],
   [ 1,  0,  0,  1,  0,  0],
   [ 0,  1,  0,  1,  0,  0],
   [ 1,  0,  0,  0,  1,  0],
   [ 0,  1,  0,  0,  1,  0],
   [ 0,  0,  1,  1,  0,  0],
   [ 1,  0,  0,  0,  0,  1],
   [ 0,  0,  1,  0,  1,  0],
   [ 0,  1,  0,  0,  0,  1],
   [ 0,  0,  1,  0,  0,  1],
 ]
 used_coeffs = []
 for i in range(len(coeffs)):
   if not (exclude_bins and i in exclude_bins):
     used_coeffs.append(coeffs[i])
 return np.array(used_coeffs)

def calculate_uncertainties(M, deltab):
  uncs = []
  for i in range(6):
    uncs.append(0)
    for k in range(len(M[0])):
      uncs[i] += (M[i, k] * deltab[k]) ** 2
    uncs[i] = math.sqrt(uncs[i])
  return np.array(uncs)

#Calculates the results
def calculate(catRatios, exclude_bins = None, weighted = True):
  A = make_coefficient_matrix(exclude_bins)
  (b, W) = make_category_matrix(catRatios, weighted)
  w = np.diag(1/W**2)
  Aw = np.dot(w, A)
  M = calculate_M(Aw)
  bw = np.dot(w, b)
  rates = calculate_rates(M, bw)
  uncs = calculate_uncertainties(np.dot(M,w), W)
  return (rates.tolist(), uncs.tolist())

def get_latex_line():
  return '\n' + '%' * 120 + '\n'

def get_solution_latex(rates, uncs, datastring):
  latex = """
\multirow{2}{*}{%s}   & $0     \leq\eta < 1.479$ & $%.4f \pm %.4f$ & $%.4f \pm %.4f$ & $%.4f \pm %.4f$  \\\\
                %s    & $1.479 \leq\eta < 2.5$   & $%.4f \pm %.4f$ & $%.4f \pm %.4f$ & $%.4f \pm %.4f$  \\\\
\hline""" % (datastring,
    rates[0] * 100, uncs[0] * 100,
    rates[1] * 100, uncs[1] * 100,
    rates[2] * 100, uncs[2] * 100,
    ' ' * len(datastring),
    rates[3] * 100, uncs[3] * 100,
    rates[4] * 100, uncs[4] * 100,
    rates[5] * 100, uncs[5] * 100,
   )
  return get_latex_line() + latex + get_latex_line()

def get_ratios_latex(ratios, datastring):
  latex = """
\multirow{2}{*}{%s} & $0     \leq \eta < 1.479$ & $%.4f \pm %.4f$ & $%.4f \pm %.4f$ & $%.4f \pm %.4f$  \\\\
                %s  & $1.479 \leq \eta < 2.5$   & $%.4f \pm %.4f$ & $%.4f \pm %.4f$ & $%.4f \pm %.4f$  \\\\
\hline""" % (datastring,
    ratios[0][0] * 100, ratios[0][1] * 100,
    ratios[1][0] * 100, ratios[1][1] * 100,
    ratios[2][0] * 100, ratios[2][1] * 100,
    ' ' * len(datastring),
    ratios[3][0] * 100, ratios[3][1] * 100,
    ratios[4][0] * 100, ratios[4][1] * 100,
    ratios[5][0] * 100, ratios[5][1] * 100,
   )
  return get_latex_line() + latex + get_latex_line()

def calculate_solution(categoryRatios, exclude_bins, output_filename, fallback_value, datastring = ""):
  rates, uncs = calculate(categoryRatios, exclude_bins)
  fit_results_to_file(rates, uncs, output_filename, fallback_value)
  if datastring:
    print(get_solution_latex(rates, uncs, datastring))
