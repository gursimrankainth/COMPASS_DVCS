import ROOT
import numpy as np
import array as arr
import pickle
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from dvcs_constants_2016 import * 
from root_likelihood import * 

# **********************************
# returns a dictionary; first index is period, second is value (sigma or err), third is t bin 
# ex. results["P04"]["sigma_muPlus"][0] -> mu+ cross section for P04 in the first t bin 
def load_cross_sections(filename):
  results = {}
  with open(filename, "r") as f:
    lines = f.readlines()
  
  i = 0
  while i < len(lines):
    line = lines[i].strip()
    if line == "":
      i += 1
      continue

    period = line.split()[0]
    results[period] = {}

    results[period]["sigma_muPlus"] = np.array(eval(line.split(":")[1]))

    i += 1
    line = lines[i].strip()
    results[period]["sigma_muMinus"] = np.array(eval(line.split(":")[1]))

    i += 1
    line = lines[i].strip()
    results[period]["err_muPlus"] = np.array(eval(line.split(":")[1]))

    i += 1
    line = lines[i].strip()
    results[period]["err_muMinus"] = np.array(eval(line.split(":")[1]))

    i += 1

  return results


# **********************************
# Get the cross section dictionary; first index is period, second is value (sigma or err), third is t bin 
# ex. results["P04"]["sigma_muPlus"][0] -> mu+ cross section for P04 in the first t bin 
with open("dvcs_xSection_results_515.pkl", "rb") as f:
  dvcs_results = pickle.load(f)

# **********************************
# Average across all data (per t-bin)
def average_charge(results, charge="muPlus"):
  if charge == "muPlus": 
    key_sigma = "sigma_muPlus"
    key_err   = "err_muPlus"
  elif charge == "muMinus":
    key_sigma = "sigma_muMinus"
    key_err   = "err_muMinus"
  else: 
    raise ValueError("Invalid charge. Use muPlus or muMinus")

  n_t = 4  # number of t bins 

  sigma_avg = []
  err_avg = []

  for t in range(n_t):
    values = []
    errors = []

    for p in PERIODS:
      values.append(results[p][key_sigma][t])
      errors.append(results[p][key_err][t])

    values = np.array(values)
    errors = np.array(errors)

    # simple (unweighted) average
    avg = np.mean(values)

    # propagate error for unweighted mean
    # assumes uncorrelated uncertainties
    err = np.sqrt(np.sum(errors**2)) / len(errors)

    sigma_avg.append(avg)
    err_avg.append(err)

  return np.array(sigma_avg), np.array(err_avg)

# **********************************
# Average across both beam charges
def average_all(sigma1, sigma2, err1, err2):
  sigma_avg = (sigma1 + sigma2) / 2
  err_avg = (err1 + err2) / 2 
  print(sigma_avg)
  print(err_avg)
  return sigma_avg, err_avg


# **********************************
# Slope fit function - reduced chi2 linear fit without bin integration
def dvcs_slope_linear(t_edges, sigma_total, err_total, tag="P04"):
  t_edges = np.array(t_edges, dtype=np.float64)
  t_centers = 0.5 * (t_edges[1:] + t_edges[:-1])
  
  t_fit = t_centers
  sigma_fit = sigma_total
  err_fit = err_total
  
  def model(t, sigma0, B):
    return sigma0 * np.exp(-B * t)
  
  sigma0_init = sigma_fit[0] if len(sigma_fit) > 0 else 1.0
  B_init = 5.0
  
  popt, pcov = curve_fit(
    model,
    t_fit,
    sigma_fit,
    sigma=err_fit,
    p0=[sigma0_init, B_init],
    absolute_sigma=True
  )

  sigma0_fit, B_fit = popt
  sigma0_err, B_err = np.sqrt(np.diag(pcov))
  
  # Plot
  plt.figure()
  plt.errorbar(t_fit, sigma_fit, yerr=err_fit, fmt='o', capsize=3, elinewidth=1)
  t_smooth = np.linspace(min(t_fit), max(t_fit), 200)
  plt.plot(t_smooth, model(t_smooth, sigma0_fit, B_fit), color='black')

  plt.xlabel("|t| (GeV$^2$)")
  plt.ylabel("dsigma/dt (nb/GeV$^2$)")
  plt.title(f"DVCS t-slope fit ({tag})")
  plt.yscale('log')

  plt.text(
    0.95, 0.95,
    f"B = {B_fit:.3f} ± {B_err:.3f} GeV$^{{-2}}$",
    transform=plt.gca().transAxes,
    verticalalignment='top',
    horizontalalignment='right',
    fontsize=10,
    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
  )

  plt.tight_layout()
  plt.savefig(f"dvcs_t_fit_{tag}.png", dpi=300)
  plt.close()
  
  print(f"Lin: B = {B_fit:.3f} ± {B_err:.3f} GeV^-2")
  #print(f"Saved plot: dvcs_t_fit_{tag}.png")
  
  return B_fit, B_err, sigma0_fit, sigma0_err


# Slope fit function - reduced chi2 exponential fit with bin-integration
def dvcs_slope_binned(t_edges, sigma, sigma_err, x0=(1.0, 5.0)):
  t_edges = np.array(t_edges, dtype=np.float64)
  sigma = np.array(sigma, dtype=np.float64)
  sigma_err = np.array(sigma_err, dtype=np.float64)

  # Bin widths
  dt = t_edges[1:] - t_edges[:-1]

  # Convert measured differential cross section to bin-integrated quantity
  sigma_int = sigma * dt
  sigma_err_int = sigma_err * dt

  # Bin edges
  t_min = t_edges[:-1]
  t_max = t_edges[1:]

  # Bin-integrated exponential model: A * exp(-B t)
  def model(params):
    A, B = params
    if np.abs(B) < 1e-12:
      return A * (t_max - t_min)
    return A * (np.exp(-B * t_min) - np.exp(-B * t_max)) / B

  # Chi-square in integrated space
  def chi2(params):
    A, B = params
    if B < 0:
      return 1e30  # enforce physical slope
    pred = model(params)
    return np.sum((sigma_int - pred)**2 / sigma_err_int**2)

  # Fit
  res = minimize(
    chi2,
    x0=np.array(x0, dtype=np.float64),
    method="BFGS"
  )

  A_fit, B_fit = res.x

  # Error estimate (Hessian approximation from BFGS)
  try:
    cov = res.hess_inv
    B_err = np.sqrt(cov[1, 1])
  except Exception:
    B_err = np.nan

  # Final chi2 / ndf
  pred = model(res.x)
  chi2_val = np.sum((sigma_int - pred)**2 / sigma_err_int**2)
  ndf = len(sigma) - 2

  print(f"Exp: B = {B_fit:.3f} ± {B_err:.3f}")

  return B_fit, B_err


# Slope fit function - binned likelihood fit using SciPy
# This function follows Johannes' procedure 
def dvcs_slope_likelihood_python(t_edges, dsigma_dt, dsigma_dt_err):
  # Bin-integrated cross section values
  N = np.zeros(len(dsigma_dt))
  N_err = np.zeros(len(dsigma_dt))

  for i in range(len(dsigma_dt)):
    delta_t = t_edges[i+1] - t_edges[i]
    N[i] = dsigma_dt[i] * delta_t
    N_err[i] = dsigma_dt_err[i] * delta_t**2

  # Exponential model per bin 
  def model_per_bin(B, t_edges):
    t_min = t_edges[0]
    t_max = t_edges[-1]

    denom = np.exp(-B * t_min) - np.exp(-B * t_max)
    mu = np.zeros(len(t_edges) - 1)

    for i in range(len(t_edges) - 1):
      num = np.exp(-B * t_edges[i]) - np.exp(-B * t_edges[i+1])
      mu[i] = num / denom

    return mu

  # Exponential likelihood 
  def neg_log_likelihood(B_array):
    B = B_array[0]
    mu = model_per_bin(B, t_edges)
    mu = np.clip(mu, 1e-12, None)
    return -np.sum(N * np.log(mu))

  # Variance likelihood
  def neg_log_likelihood_var(B_array):
    B = B_array[0]
    sigma_t_error = np.zeros(len(dsigma_dt_err))

    for i in range(len(dsigma_dt_err)):
      delta_t = t_edges[i+1] - t_edges[i]
      sigma_t_error[i] = dsigma_dt_err[i] * delta_t**2

    t_min = t_edges[0]
    t_max = t_edges[-1]

    denom = np.exp(-B * t_min) - np.exp(-B * t_max)
    norm_factor = np.sum(sigma_t_error) / denom

    mu = np.zeros(len(t_edges) - 1)
    for i in range(len(t_edges) - 1):
      mu[i] = norm_factor * (
        np.exp(-B * t_edges[i]) - np.exp(-B * t_edges[i+1])
      )

    mu = np.clip(mu, 1e-12, None)
    return -np.sum(sigma_t_error * np.log(mu))

  # Fit main expoeneital function 
  res = minimize(
    neg_log_likelihood,
    x0=np.array([1.0]),
    method="BFGS"
  )

  B_fit = res.x[0]
  cov_main = res.hess_inv[0, 0]

  # Fit variance function 
  res_var = minimize(
    neg_log_likelihood_var,
    x0=np.array([1.0]),
    method="BFGS"
  )

  B_var_fit = res_var.x[0]
  cov_var = res_var.hess_inv[0, 0]

  # Combined error (this is the final error that is reported)
  B_err_combined = np.sqrt(cov_main**2 / cov_var)

  print(f"LLP: B = {B_fit:.3f} ± {B_err_combined:.3f} GeV^-2")
  return B_fit, B_err_combined

# **********************************
# Average slope using per period calculated values 
def weighted_avg_slope(values, errors):
  weights = 1 / errors**2
  avg = np.sum(values * weights) / np.sum(weights)
  err = 1 / np.sqrt(np.sum(weights))
  return avg, err

# **********************************
# Main
#results = load_cross_sections("dvcs_xSecVals.txt")
#print(dvcs_results["P04"]["sigma_muPlus"])
t_edges = [0.08, 0.136, 0.219, 0.36, 0.64]
 
# ****************
# Average across all periods
sigma_muPlus_avg, err_muPlus_avg = average_charge(dvcs_results, charge="muPlus")
sigma_muMinus_avg, err_muMinus_avg = average_charge(dvcs_results, charge="muMinus")
print("\n--- Avg. cross sections per beam charge ---")
print("mu+ dsigma/dt (nb/GeV²):", sigma_muPlus_avg)
print("mu+ error (nb/GeV²)", err_muPlus_avg)
print("mu- dsigma/dt (nb/GeV²):", sigma_muMinus_avg)
print("mu- error (nb/GeV²)", err_muMinus_avg)

# ****************
# t-slopes per beam charge 
print("\n--- Slopes per beam charge ---")
dvcs_slope_likelihood_root(t_edges, sigma_muPlus_avg, err_muPlus_avg)
dvcs_slope_likelihood_root(t_edges, sigma_muMinus_avg, err_muMinus_avg)

# ****************
sigma_combined, err_combined = average_all(sigma_muPlus_avg, sigma_muMinus_avg, err_muPlus_avg, err_muMinus_avg)
print("\n--- Avg. cross section ---")
print("dsigma/dt (nb/GeV²):", sigma_combined)
print("error (nb/GeV²)", err_combined)

print("\n--- Slope ---")
B_fit, B_var, B_err = dvcs_slope_likelihood_root(t_edges, sigma_combined, err_combined)
_, _, sigma0_fit, sigma0_err = dvcs_slope_linear(t_edges, sigma_combined, err_combined)

# ****************
# Per period average cross section 
P04, P04_err = average_all(dvcs_results["P04"]["sigma_muPlus"], dvcs_results["P04"]["sigma_muMinus"], dvcs_results["P04"]["err_muPlus"], dvcs_results["P04"]["err_muMinus"])
P05, P05_err = average_all(dvcs_results["P05"]["sigma_muPlus"], dvcs_results["P05"]["sigma_muMinus"], dvcs_results["P05"]["err_muPlus"], dvcs_results["P05"]["err_muMinus"])
P06, P06_err = average_all(dvcs_results["P06"]["sigma_muPlus"], dvcs_results["P06"]["sigma_muMinus"], dvcs_results["P06"]["err_muPlus"], dvcs_results["P06"]["err_muMinus"])
P07, P07_err = average_all(dvcs_results["P07"]["sigma_muPlus"], dvcs_results["P07"]["sigma_muMinus"], dvcs_results["P07"]["err_muPlus"], dvcs_results["P07"]["err_muMinus"])
P08, P08_err = average_all(dvcs_results["P08"]["sigma_muPlus"], dvcs_results["P08"]["sigma_muMinus"], dvcs_results["P08"]["err_muPlus"], dvcs_results["P08"]["err_muMinus"])
P09, P09_err = average_all(dvcs_results["P09"]["sigma_muPlus"], dvcs_results["P09"]["sigma_muMinus"], dvcs_results["P09"]["err_muPlus"], dvcs_results["P09"]["err_muMinus"])
""" """
# ****************
# Test with Johannes' values 
""" print("\n--- Test with 2023 results ---")
J_muPlus_avg = np.array([24.646, 16.894, 10.403, 3.655])
J_muPlus_err_avg = np.array([4.092, 2.505, 1.505, 0.737])
J_muMinus_avg = np.array([31.907, 18.010, 9.691, 2.54])
J_muMinus_err_avg = np.array([4.538, 2.930, 1.677, 0.744])
#dvcs_slope_likelihood_python(t_edges, J_muPlus_avg, J_muPlus_err_avg)
#dvcs_slope_likelihood_python(t_edges, J_muMinus_avg, J_muMinus_err_avg) 

J_combined, J_err_combined = average_all(J_muPlus_avg, J_muMinus_avg, J_muPlus_err_avg, J_muMinus_err_avg)
_, _, sigma0_fit, sigma0_err = dvcs_slope_linear(t_edges, J_combined, J_err_combined)
J_B_fit, J_B_var_fit, J_B_err_combined = dvcs_slope_likelihood_python(t_edges, J_combined, J_err_combined) """

# ****************
# Make the plot 
# Prep the data for the plot  
data = [
  (P04, P04_err, 'P04'),
  (P05, P05_err, 'P05'),
  (P06, P06_err, 'P06'),
  (P07, P07_err, 'P07'),
  (P08, P08_err, 'P08'),
  (P09, P09_err, 'P09')
] 
""" """ 

def make_plot(data, B, B_err, y_int):
  # t bin centers
  t_cent = np.array([0.10, 0.18, 0.28, 0.48])

  # offsets for visual separation
  offsets = [-0.01, -0.006, -0.002, -0.004, -0.01, 0.0, 0.0025]

  plt.figure()

  # plot ONLY once (offset version)
  for (y, err, label), dx in zip(data, offsets):
    plt.errorbar(t_cent + dx, y, yerr=err, fmt='o', capsize=4, label=label)

  t_line = np.linspace(t_cent.min(), t_cent.max(), 200)
  y_line = y_int * np.exp(-B * t_line)
  y_up = y_int * np.exp(-(B - B_err) * t_line)
  y_down = y_int * np.exp(-(B + B_err) * t_line)

  plt.plot(t_line, y_line, 'k-', label=rf'$e^{{-Bt}},\ B = {B:.3f} \pm {B_err:.3f}$')
  plt.fill_between(t_line, y_down, y_up, color='green', alpha=0.25, label='B uncertainty')

  plt.xlabel(r'$|t|\;(\mathrm{GeV}/c)^2$')
  plt.ylabel(r'$\frac{d\sigma}{d|t|}\;[\mathrm{nb}\,(\mathrm{GeV}/c)^{-2}]$')
  plt.yscale('log') 
  plt.legend()

  plt.savefig("test_plot.png", dpi=300, bbox_inches="tight")

make_plot(data=data, B=B_fit, B_err=B_err, y_int=sigma0_fit)
