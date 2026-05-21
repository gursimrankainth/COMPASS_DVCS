import ROOT
from ROOT import TLorentzVector
import pickle
import os 
import math
import numpy as np
import matplotlib.pyplot as plt
from dvcs_constants_2016 import *
from scipy.optimize import curve_fit

# *******************************************************************
# *                    *** DATA PREP ***                            *
# *******************************************************************
# **********************************
# Constants 
alpha_em = 0.0072973525693 # electromagnetic fine structure constant 
M_mu = 105.6583755e-3 # GeV/c

# **********************************
# Real data 
real_dir = "/eos/user/g/gkainth/real/"
real_files = [os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root")]

# HEPGEN BH MC Data (Generated data)
gen_dir = "/eos/user/g/gkainth/BH/"
gen_files = [os.path.join(gen_dir, "gen_P04_muPlus.root"), os.path.join(gen_dir, "gen_P04_muMinus.root"),
             os.path.join(gen_dir, "gen_P05_muPlus.root"), os.path.join(gen_dir, "gen_P05_muMinus.root"),
             os.path.join(gen_dir, "gen_P06_muPlus.root"), os.path.join(gen_dir, "gen_P06_muMinus.root"),
             os.path.join(gen_dir, "gen_P07_muPlus.root"), os.path.join(gen_dir, "gen_P07_muMinus.root"),
             os.path.join(gen_dir, "gen_P08_muPlus.root"), os.path.join(gen_dir, "gen_P08_muMinus.root"),
             os.path.join(gen_dir, "gen_P09_muPlus.root"), os.path.join(gen_dir, "gen_P09_muMinus.root")]

# HEPGEN BH MC Data (Reconstructed data)
hepBH_dir = "/eos/user/g/gkainth/BH/"
hepBH_files = [os.path.join(hepBH_dir, "filtered_P04_muPlus.root"), os.path.join(hepBH_dir, "filtered_P04_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P05_muPlus.root"), os.path.join(hepBH_dir, "filtered_P05_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P06_muPlus.root"), os.path.join(hepBH_dir, "filtered_P06_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P07_muPlus.root"), os.path.join(hepBH_dir, "filtered_P07_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P08_muPlus.root"), os.path.join(hepBH_dir, "filtered_P08_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P09_muPlus.root"), os.path.join(hepBH_dir, "filtered_P09_muMinus.root")]

# HEPGEN Invisible Pi0 MC Data (Reconstructed data)
hepPi0_dir = "/eos/user/g/gkainth/HepgenPi0/"
hepPi0_files = [os.path.join(hepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P09_muMinus.root")]

# LEPTO Invisible Pi0 MC Data (Reconstructed data)
lepPi0_dir = "/eos/user/g/gkainth/LeptoPi0/"
lepPi0_files = [os.path.join(lepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P09_muMinus.root")]

# **********************************
# Total Luminosity 
tot_lum_muPlus = np.sum(LUMINOSITY_MUPLUS)
tot_lum_muMinus = np.sum(LUMINOSITY_MUMINUS)

# **********************************
# Binning scheme (same as what is used for the final accpetance)
# nu: 4 bins of width 5.5 GeV between 10 and 32 GeV
nu_edges = np.linspace(10, 32, 5)
nu_bins = list(zip(nu_edges[:-1], nu_edges[1:]))
# Q2: 4 bins of width 1 (GeV/c)^2 between 1 and 5
Q2_edges = np.linspace(1, 5, 5)
Q2_bins = list(zip(Q2_edges[:-1], Q2_edges[1:]))
# |t|: 4 bins -> each bin should have roughly the same no. of events 
t_edges = [0.08, 0.136, 0.219, 0.36, 0.64]
t_bins = list(zip(t_edges[:-1], t_edges[1:]))
# phi: 8 bins of width pi/4 rad between -pi and pi
phi_edges = np.linspace(-np.pi, np.pi, 9)
phi_bins = list(zip(phi_edges[:-1], phi_edges[1:]))

# Compute the bin widths
# bin size is uneven in t so it cannot be calculated like this 
delta_nu   = np.diff(nu_edges)       # array of size n_nu_bins
delta_Q2   = np.diff(Q2_edges)       # array of size n_Q2_bins
delta_phi  = np.diff(phi_edges)      # array of size n_phi_bins


# *******************************************************************
# *                    *** ACCEPTANCE ***                           *
# *******************************************************************
# **********************************
# Fill arrays with sum of weights in each bin 
def fill_weights(data_type="gen", period="P04"):
  # Fixed binning structure (Axis order: [t][Q2][nu][phi])
  shape = (4, 4, 4, 8)
  weights_muPlus     = np.zeros(shape, dtype=np.float64)
  weights_muMinus    = np.zeros(shape, dtype=np.float64)
  weights_muPlus_sq  = np.zeros(shape, dtype=np.float64)
  weights_muMinus_sq = np.zeros(shape, dtype=np.float64)

  n_t, n_Q2, n_nu, n_phi = shape

  # Select tree and files
  if data_type == "gen":
    chain = ROOT.TChain("USR970_gen")
    file_list = gen_files

  elif data_type == "rec":
    chain = ROOT.TChain("USR970_filtered")
    file_list = hepBH_files

  else:
    raise ValueError("data_type must be 'gen' or 'rec'")

  for f in file_list:
    if period in f:
      chain.Add(f)

  # Event loop
  for event in chain:
    charge = event.Q_beam
    weight = event.weight_DVCS

    if data_type == "gen":
      Q2  = event.Q2_gen
      nu  = event.nu_gen
      t   = event.t_gen
      phi = event.phi_gg_gen
    else:
      Q2  = event.Q2_fit
      nu  = event.nu_fit
      t   = event.t_fit
      phi = event.phi_gg_fit

    # Normalize phi to [-pi, pi)
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    # Bin Q2
    Q2_bin  = np.searchsorted(Q2_edges, Q2) - 1
    if Q2_bin < 0 or Q2_bin >= n_Q2:
      continue

    # Bin nu
    nu_bin  = np.searchsorted(nu_edges, nu) - 1
    if nu_bin < 0 or nu_bin >= n_nu:
      continue

    # Bin phi
    phi_bin = np.searchsorted(phi_edges, phi) - 1
    if phi_bin < 0 or phi_bin >= n_phi:
      continue

    # Bin |t|
    t_abs = abs(t)
    t_bin = -1
    for i, (t_low, t_high) in enumerate(t_bins):
      if t_low <= t_abs < t_high:
        t_bin = i
        break

    if t_bin == -1:
      continue

    # Fill arrays
    if charge == 1:
      weights_muPlus[t_bin, Q2_bin, nu_bin, phi_bin]     += weight
      weights_muPlus_sq[t_bin, Q2_bin, nu_bin, phi_bin]  += weight * weight

    elif charge == -1:
      weights_muMinus[t_bin, Q2_bin, nu_bin, phi_bin]    += weight
      weights_muMinus_sq[t_bin, Q2_bin, nu_bin, phi_bin] += weight * weight

  return (weights_muPlus, weights_muMinus, weights_muPlus_sq, weights_muMinus_sq)

# **********************************
# Get the acceptance 
def compute_acceptance(rec_muPlus, rec_muMinus, gen_muPlus, gen_muMinus):
  shape = rec_muPlus.shape

  acc_muPlus  = np.zeros(shape, dtype=np.float64)
  acc_muMinus = np.zeros(shape, dtype=np.float64)

  # muPlus 
  mask_plus = gen_muPlus != 0.0
  acc_muPlus[mask_plus] = (
    rec_muPlus[mask_plus] / gen_muPlus[mask_plus]
  )

  # muMinus
  mask_minus = gen_muMinus != 0.0
  acc_muMinus[mask_minus] = (
    rec_muMinus[mask_minus] / gen_muMinus[mask_minus]
  )

  return acc_muPlus, acc_muMinus

# **********************************
# Get the acceptance error (returns variance and standard deviation)
def compute_acceptance_error(rec_muPlus_sq, rec_muMinus_sq,
                             gen_muPlus_sq, gen_muMinus_sq,
                             rec_muPlus, rec_muMinus,
                             gen_muPlus, gen_muMinus):

  shape = rec_muPlus.shape

  var_muPlus  = np.zeros(shape, dtype=np.float64)
  var_muMinus = np.zeros(shape, dtype=np.float64)

  err_muPlus  = np.zeros(shape, dtype=np.float64)
  err_muMinus = np.zeros(shape, dtype=np.float64)

  for idx in np.ndindex(shape):
    # μ+
    weightRec   = rec_muPlus[idx]
    weightGen   = gen_muPlus[idx]
    weight2Rec  = rec_muPlus_sq[idx]
    weight2Gen  = gen_muPlus_sq[idx]

    if weightGen == 0.0:
      varAcc = 0.0
    else:
      errT1 = (1.0 / weightGen)**2 * weight2Rec
      errT2 = (weightRec / (weightGen**2))**2 * weight2Gen
      varAcc = errT1 + errT2

    var_muPlus[idx] = varAcc
    err_muPlus[idx] = np.sqrt(varAcc)

    # μ-
    weightRec   = rec_muMinus[idx]
    weightGen   = gen_muMinus[idx]
    weight2Rec  = rec_muMinus_sq[idx]
    weight2Gen  = gen_muMinus_sq[idx]

    if weightGen == 0.0:
      varAcc = 0.0
    else:
      errT1 = (1.0 / weightGen)**2 * weight2Rec
      errT2 = (weightRec / (weightGen**2))**2 * weight2Gen
      varAcc = errT1 + errT2

    var_muMinus[idx] = varAcc
    err_muMinus[idx] = np.sqrt(varAcc)

  return var_muPlus, var_muMinus, err_muPlus, err_muMinus

# **********************************
# Check that the acceptance and associated error make sense (test by integrating over phi)
# Plotting funtion is used in test_acceptance
def plot_acceptance_integrated_phi(acceptance_muPlus, acceptance_err_muPlus, 
                                   acceptance_muMinus, acceptance_err_muMinus):
  nu_bin_centers = 0.5 * (nu_edges[:-1] + nu_edges[1:])
  n_t_bins, n_Q2_bins, n_nu_bins = acceptance_muPlus.shape

  fig, axes = plt.subplots(nrows=n_Q2_bins, ncols=n_t_bins, figsize=(18, 22), sharex=True, sharey=True)

  for i in range(n_Q2_bins):  # rows (Q2 -> y values)    
    for j in range(n_t_bins):  # columns (t -> x)   
      ax = axes[i, j]

      # acceptance_muPlus and acceptance_muMinus shape: (t, Q2, nu)
      y_muPlus = acceptance_muPlus[j, i]
      yerr_muPlus = acceptance_err_muPlus[j, i]
      y_muMinus = acceptance_muMinus[j, i]
      yerr_muMinus = acceptance_err_muMinus[j, i]

      ax.errorbar(nu_bin_centers, y_muPlus, yerr=yerr_muPlus, fmt='o', markerfacecolor='none', 
                  markeredgecolor='red', markersize=5, ecolor='red', label='μ⁺' if i == 0 and j == 0 else "")
      ax.errorbar(nu_bin_centers, y_muMinus, yerr=yerr_muMinus, fmt='o', markerfacecolor='none', 
                  markeredgecolor='black', markersize=5, ecolor='black', label='μ⁻' if i == 0 and j == 0 else "")

      ax.set_ylim(0, 0.8)
      ax.grid(True, linestyle='--', linewidth=0.5)
      ax.axhline(0, color='gray', linewidth=0.5)

  # Legend
  handles, labels = axes[0, 0].get_legend_handles_labels()
  fig.legend(handles, labels, loc='upper right', fontsize=14, markerscale=1.5)

  # Create a new set of axes for the nu and acceptance scale at the top right corner
  nu_axis = fig.add_axes([0.71, 0.92, 0.185, 0.035])  # [left, bottom, width, height]
  acc_axis = fig.add_axes([0.8, 0.86, 0.11, 0.0825])

  # Setup nu axis
  nu_axis.set_xlim(10, 32)
  nu_axis.set_xticks(nu_edges)
  nu_axis.set_xticklabels([f"{edge:.0f}" for edge in nu_edges], fontsize=12)
  nu_axis.set_yticks([])
  nu_axis.yaxis.set_visible(False)
  nu_axis.tick_params(axis='x', direction='in', length=5, top=True, bottom=False)
  nu_axis.xaxis.set_label_position('top')
  nu_axis.set_xlabel(r"$\nu$ [GeV]", fontsize=14, labelpad=20)
  nu_axis.xaxis.tick_top()
  nu_axis.patch.set_facecolor('none')
  for name, spine in nu_axis.spines.items():
    spine.set_visible(name == 'top')
    if name == 'top':
      spine.set_linewidth(1.0)
      spine.set_color('black')

  # Setup acceptance axis
  acc_axis.set_ylim(0, 0.8)
  acc_axis.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
  acc_axis.set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8"], fontsize=12)
  acc_axis.set_xticks([])
  acc_axis.xaxis.set_visible(False)
  acc_axis.set_ylabel("Acceptance", fontsize=14, labelpad=15)
  acc_axis.yaxis.set_label_position('right')
  acc_axis.yaxis.tick_right()
  acc_axis.tick_params(axis='y', direction='in', length=5)
  acc_axis.patch.set_facecolor('none')
  for name, spine in acc_axis.spines.items():
    spine.set_visible(name == 'right')
    if name == 'right':
      spine.set_linewidth(1.0)
      spine.set_color('black')

  # Global axes for t and Q2 
  t_axis = fig.add_axes([0.125, 0.09, 0.775, 0.035])
  Q2_axis = fig.add_axes([0.083, 0.11, 0.034, 0.835])

  # t axis setup
  tick_positions = np.arange(len(t_edges))
  t_axis.set_xlim(0, len(t_edges) - 1)
  t_axis.set_xticks(tick_positions)
  t_axis.set_xticklabels([f"{a}" for a in t_edges], fontsize=12)
  t_axis.set_yticks([])
  t_axis.xaxis.tick_bottom()
  t_axis.set_xlabel(r"$|t|$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
  t_axis.patch.set_facecolor('none')
  for name, spine in t_axis.spines.items():
    spine.set_visible(name == 'bottom')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  # Q2 axis setup
  Q2_axis.set_ylim(Q2_edges[-1], Q2_edges[0])
  Q2_axis.set_yticks(Q2_edges)
  Q2_axis.set_yticklabels([f"{int(a)}" for a in Q2_edges], fontsize=12)
  Q2_axis.set_xticks([])
  Q2_axis.yaxis.tick_left()
  Q2_axis.set_ylabel(r"$Q^2$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
  Q2_axis.patch.set_facecolor('none')
  for name, spine in Q2_axis.spines.items():
    spine.set_visible(name == 'left')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  plt.tight_layout(rect=[0.10, 0.10, 0.9, 0.95])
  plt.savefig("acceptance_integrated_phi.png", dpi=300)

# Integrate over phi and plot the acceptance 
def test_acceptance(rec_muPlus, rec_muMinus, rec_muPlus_sq, rec_muMinus_sq,
                    gen_muPlus, gen_muMinus, gen_muPlus_sq, gen_muMinus_sq):
    
    # Sum over phi (same logic as in acceptance.py)
    rec_muPlus_phiInt  = np.sum(rec_muPlus, axis=3)
    rec_muMinus_phiInt = np.sum(rec_muMinus, axis=3)
    gen_muPlus_phiInt  = np.sum(gen_muPlus, axis=3)
    gen_muMinus_phiInt = np.sum(gen_muMinus, axis=3)

    rec_muPlus_sq_phiInt  = np.sum(rec_muPlus_sq, axis=3)
    rec_muMinus_sq_phiInt = np.sum(rec_muMinus_sq, axis=3)
    gen_muPlus_sq_phiInt  = np.sum(gen_muPlus_sq, axis=3)
    gen_muMinus_sq_phiInt = np.sum(gen_muMinus_sq, axis=3)
    print("Finished integrating over phi for all arrays")

    # Acceptance
    acceptance_muPlus = np.zeros_like(gen_muPlus_phiInt)
    acceptance_muMinus = np.zeros_like(gen_muMinus_phiInt)

    nonzero_muPlus = gen_muPlus_phiInt != 0
    nonzero_muMinus = gen_muMinus_phiInt != 0

    acceptance_muPlus[nonzero_muPlus] = rec_muPlus_phiInt[nonzero_muPlus] / gen_muPlus_phiInt[nonzero_muPlus]
    acceptance_muMinus[nonzero_muMinus] = rec_muMinus_phiInt[nonzero_muMinus] / gen_muMinus_phiInt[nonzero_muMinus]
    print("Finished computing acceptance")

    # Error
    var_muPlus, var_muMinus, err_muPlus, err_muMinus = compute_acceptance_error(
      rec_muPlus_sq_phiInt, rec_muMinus_sq_phiInt,
      gen_muPlus_sq_phiInt, gen_muMinus_sq_phiInt,
      rec_muPlus_phiInt, rec_muMinus_phiInt,
      gen_muPlus_phiInt, gen_muMinus_phiInt
    )
    print("Finished computing acceptance error")

    # Plot
    plot_acceptance_integrated_phi(acceptance_muPlus, err_muPlus, acceptance_muMinus, err_muMinus)
    print("Acceptance plot created successfully!")


# *******************************************************************
# *                    *** CROSS SECTION ***                        *
# *******************************************************************
# **********************************
# Transverse virtual photon flux factor 
def getfluxFac(Q2, nu, y, E_mu, xbj):
  c1 = (alpha_em / (2 * np.pi)) * ((1 - xbj) / (Q2 * y * E_mu))
  c2 = y**2 * (1 - (2 * M_mu**2 / Q2))
  c3 = (2 / (1 + Q2 / nu**2)) * (1 - y - Q2 / (4 * E_mu**2))
  virFluxFac = c1 * (c2 + c3)
  return virFluxFac

# **********************************
# Binned sum over the unweighted data for a single period (use for real and LEPTO pi0)
def unweighted_sum(muPlus_array, muMinus_array, files, data, period="P04"): 
  if data == "real":
    chain = ROOT.TChain("USR970_filtered")
  elif data == "lepPi0":
    chain = ROOT.TChain("USR970_filtered")
  else: 
    raise ValueError("data must be 'real' or 'lepPi0'")

  for f in files: 
    if period in f:
      chain.Add(f)

  for event in chain: 
    charge = event.Q_beam
    Q2 = event.Q2_fit
    nu = event.nu_fit
    y = event.y_fit
    t = event.t_fit
    xbj = event.xbj_fit
    phi = event.phi_gg_fit
    E_inMu = event.inMuFit_TL.E()

    # Normalize phi 
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    # Get the virtual photon flux factor
    fluxFac = getfluxFac(Q2, nu, y, E_inMu, xbj)

    # Bin in Q2
    Q2_bin  = np.searchsorted(Q2_edges, Q2) - 1
    if Q2_bin < 0 or Q2_bin >= len(Q2_bins): 
      continue

    # Bin nu
    nu_bin  = np.searchsorted(nu_edges, nu) - 1
    if nu_bin < 0 or nu_bin >= len(nu_bins):
      continue

    # Bin phi
    phi_bin = np.searchsorted(phi_edges, phi) - 1
    if phi_bin < 0 or phi_bin >= len(phi_bins):
      continue

    # Bin |t|
    t_bin = -1
    for i, (t_low, t_high) in enumerate(t_bins):
      if t_low <= np.abs(t) < t_high:
        t_bin = i
    if t_bin == -1: 
      continue 

    # Axis order: [t][Q2][nu][phi]
    i1, i2, i3, i4 = t_bin, Q2_bin, nu_bin, phi_bin

    if charge == 1:
      muPlus_array[i1][i2][i3][i4] += (1/fluxFac)
    elif charge == -1: 
      muMinus_array[i1][i2][i3][i4] += (1/fluxFac)

# **********************************
# Binned sum over the weighted data for a single period (use for HEPGEN BH and HEPGEN pi0)
def weighted_sum(muPlus_array, muMinus_array, files, data, period="P04"): 
  if data == "hepBH":
    chain = ROOT.TChain("USR970_filtered")
  elif data == "hepPi0":
    chain = ROOT.TChain("USR970_filtered")
  else: 
    raise ValueError("data must be 'hepBH' or 'hepPi0'")

  for f in files: 
    if period in f: 
      chain.Add(f)

  for event in chain: 
    charge = event.Q_beam
    Q2 = event.Q2_fit
    nu = event.nu_fit
    y = event.y_fit
    t = event.t_fit
    xbj = event.xbj_fit
    phi = event.phi_gg_fit
    E_inMu = event.inMuFit_TL.E()
    if data == "hepBH":
      weight = event.weight_PAMBH 
    else: 
      weight = event.weight_all  

    # Normalize phi 
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    # Get the virtual photon flux factor
    fluxFac = getfluxFac(Q2, nu, y, E_inMu, xbj)

    # Bin in Q2
    Q2_bin  = np.searchsorted(Q2_edges, Q2) - 1
    if Q2_bin < 0 or Q2_bin >= len(Q2_bins): 
      continue

    # Bin nu
    nu_bin  = np.searchsorted(nu_edges, nu) - 1
    if nu_bin < 0 or nu_bin >= len(nu_bins):
      continue

    # Bin phi
    phi_bin = np.searchsorted(phi_edges, phi) - 1
    if phi_bin < 0 or phi_bin >= len(phi_bins):
      continue

    # Bin |t|
    t_bin = -1
    for i, (t_low, t_high) in enumerate(t_bins):
      if t_low <= np.abs(t) < t_high:
        t_bin = i
    if t_bin == -1: 
      continue 

    # Axis order: [t][Q2][nu][phi]
    i1, i2, i3, i4 = t_bin, Q2_bin, nu_bin, phi_bin

    if charge == 1:
      muPlus_array[i1][i2][i3][i4] += (weight/fluxFac)
    elif charge == -1: 
      muMinus_array[i1][i2][i3][i4] += (weight/fluxFac)

# **********************************
# Get sum term in 4D cross section (per period)
def get_S(real_sum, BH_sum, lepPi0_sum, hepPi0_sum, period_idx, charge="muPlus"):
  if charge == "muPlus":
    sum_term = real_sum - CBH_MUPLUS[period_idx]*BH_sum - CPI0_LEP_MUPLUS[period_idx]*R_LEPTO*lepPi0_sum - CPI0_HEP_MUPLUS[period_idx]*(1-R_LEPTO)*hepPi0_sum
  elif charge == "muMinus": 
    sum_term = real_sum - CBH_MUMINUS[period_idx]*BH_sum - CPI0_LEP_MUMINUS[period_idx]*R_LEPTO*lepPi0_sum - CPI0_HEP_MUMINUS[period_idx]*(1-R_LEPTO)*hepPi0_sum
  else:
    print('Invalid charge, please use "muPlus" or "muMinus"')
  return sum_term

# **********************************
# Get the mean cross section in (|t|, phi)
def compute_sigma_t_phi(Ncorr_ijkl, luminosity, t_edges, delta_phi, Q2_edges, nu_edges):
  n_t, n_Q2, n_nu, n_phi = Ncorr_ijkl.shape
  sigma_t_phi = np.zeros((n_t, n_phi), dtype=np.float64)

  delta_Q2_total = Q2_edges[-1] - Q2_edges[0]
  delta_nu_total = nu_edges[-1] - nu_edges[0]

  for it in range(n_t):
    delta_t = t_edges[it+1] - t_edges[it]

    for iphi in range(n_phi):
      # Sum over Q2 and nu
      S = 0.0
      for iq in range(n_Q2):
        for inu in range(n_nu):
          S += Ncorr_ijkl[it, iq, inu, iphi]

      sigma_t_phi[it, iphi] = (
        S /
        (luminosity *
         delta_t *
         delta_Q2_total *
         delta_nu_total *
         delta_phi[iphi])
      )

  return sigma_t_phi

# Integrate over phi (discrete sum over phi bins)
def integrate_over_phi(sigma_t_phi, delta_phi):
  n_t, n_phi = sigma_t_phi.shape
  sigma_t = np.zeros(n_t, dtype=np.float64)

  for it in range(n_t):
    total = 0.0
    for iphi in range(n_phi):
      total += sigma_t_phi[it, iphi] * delta_phi[iphi]
    sigma_t[it] = total

  return sigma_t


# *******************************************************************
# *                *** CROSS SECTION ERROR ***                      *
# *******************************************************************
# **********************************
# Get the (ΔD_ijkl)^2 or (ΔL_ijkl)^2 array per period 
def unweighted_err_sum(muPlus_array, muMinus_array, files, data, period="P04"): 
  # dont take the square root here, do that outside of the function
  if data == "real":
    chain = ROOT.TChain("USR970_filtered")
  elif data == "lepPi0":
    chain = ROOT.TChain("USR970_filtered")
  else: 
    raise ValueError("data must be 'real' or 'lepPi0'")

  for f in files: 
    if period in f:
      chain.Add(f)

  for event in chain: 
    charge = event.Q_beam
    Q2 = event.Q2_fit
    nu = event.nu_fit
    y = event.y_fit
    t = event.t_fit
    xbj = event.xbj_fit
    phi = event.phi_gg_fit
    E_inMu = event.inMuFit_TL.E()

    # Normalize phi 
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    # Get the virtual photon flux factor
    fluxFac = getfluxFac(Q2, nu, y, E_inMu, xbj)

    # Bin in Q2
    Q2_bin  = np.searchsorted(Q2_edges, Q2) - 1
    if Q2_bin < 0 or Q2_bin >= len(Q2_bins): 
      continue

    # Bin nu
    nu_bin  = np.searchsorted(nu_edges, nu) - 1
    if nu_bin < 0 or nu_bin >= len(nu_bins):
      continue

    # Bin phi
    phi_bin = np.searchsorted(phi_edges, phi) - 1
    if phi_bin < 0 or phi_bin >= len(phi_bins):
      continue

    # Bin |t|
    t_bin = -1
    for i, (t_low, t_high) in enumerate(t_bins):
      if t_low <= np.abs(t) < t_high:
        t_bin = i
    if t_bin == -1: 
      continue 

    # Axis order: [t][Q2][nu][phi]
    i1, i2, i3, i4 = t_bin, Q2_bin, nu_bin, phi_bin

    if charge == 1:
      muPlus_array[i1][i2][i3][i4] += (1/fluxFac)**2
    elif charge == -1: 
      muMinus_array[i1][i2][i3][i4] += (1/fluxFac)**2

# Get the (ΔB_ijkl)^2 or (ΔH_ijkl)^2 array per period
def weighted_err_sum(muPlus_array, muMinus_array, files, data, period="P04"): 
  # dont take the square root since we sqaure the sum when estimating the erro
  if data == "hepBH":
    chain = ROOT.TChain("USR970_filtered")
  elif data == "hepPi0":
    chain = ROOT.TChain("USR970_filtered")
  else: 
    raise ValueError("data must be 'hepBH' or 'hepPi0'")

  for f in files: 
    if period in f: 
      chain.Add(f)

  for event in chain: 
    charge = event.Q_beam
    Q2 = event.Q2_fit
    nu = event.nu_fit
    y = event.y_fit
    t = event.t_fit
    xbj = event.xbj_fit
    phi = event.phi_gg_fit
    E_inMu = event.inMuFit_TL.E()

    # Normalize phi 
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi

    if data == "hepBH":
      weight = event.weight_PAMBH
    else: 
      weight = event.weight_all

    # Get the virtual photon flux factor
    fluxFac = getfluxFac(Q2, nu, y, E_inMu, xbj)

    # Bin in Q2
    Q2_bin  = np.searchsorted(Q2_edges, Q2) - 1
    if Q2_bin < 0 or Q2_bin >= len(Q2_bins): 
      continue

    # Bin nu
    nu_bin  = np.searchsorted(nu_edges, nu) - 1
    if nu_bin < 0 or nu_bin >= len(nu_bins):
      continue

    # Bin phi
    phi_bin = np.searchsorted(phi_edges, phi) - 1
    if phi_bin < 0 or phi_bin >= len(phi_bins):
      continue

    # Bin |t|
    t_bin = -1
    for i, (t_low, t_high) in enumerate(t_bins):
      if t_low <= np.abs(t) < t_high:
        t_bin = i
    if t_bin == -1: 
      continue 

    # Axis order: [t][Q2][nu][phi]
    i1, i2, i3, i4 = t_bin, Q2_bin, nu_bin, phi_bin

    if charge == 1:
      muPlus_array[i1][i2][i3][i4] += (weight/fluxFac)**2
    elif charge == -1: 
      muMinus_array[i1][i2][i3][i4] += (weight/fluxFac)**2

# **********************************
# Get error for sum term (ΔS_ijkl)^2 in 4D cross section per period
def get_err_S2(D_array, B_array, L_array, H_array, period_idx, charge="muPlus"):
  # the arrays are just the sums, NOT the sqaure root of the sums (no need to sqaure them)
  t1 = D_array
  if charge == "muPlus":
    t2 = (CBH_MUPLUS[period_idx])**2 * B_array
    t3 = (CPI0_LEP_MUPLUS[period_idx])**2 * (R_LEPTO)**2 * L_array
    t4 = (CPI0_HEP_MUPLUS[period_idx])**2 * (1-R_LEPTO)**2 * H_array
    sum_term = t1 + t2 + t3 + t4
  elif charge == "muMinus": 
    t2 = (CBH_MUMINUS[period_idx])**2 * B_array
    t3 = (CPI0_LEP_MUMINUS[period_idx])**2 * (R_LEPTO)**2 * L_array
    t4 = (CPI0_HEP_MUMINUS[period_idx])**2 * (1-R_LEPTO)**2 * H_array
    sum_term = t1 + t2 + t3 + t4
  else:
    print('Invalid charge, please use "muPlus" or "muMinus"')
  return sum_term

# Get error for 4D cross section per period
def compute_4D_error(S_ijkl, A_ijkl, varS_ijkl, varA_ijkl):
    Ncorr = S_ijkl / A_ijkl
    err2 = np.zeros_like(Ncorr)

    # Mask bins with negligible signal
    mask = (S_ijkl != 0) & (A_ijkl != 0)

    # Standard variance propagation formula
    err2[mask] = (Ncorr[mask]**2) * (
      varS_ijkl[mask] / (S_ijkl[mask]**2) +
      varA_ijkl[mask] / (A_ijkl[mask]**2)
    )

    return err2

# **********************************
# Get the error for the mean cross section in |t|, phi 
def get_err_t_phi(var4D_array, luminosity, t_edges, delta_phi, Q2_edges, nu_edges):
  n_t, n_Q2, n_nu, n_phi = var4D_array.shape
  var_t_phi = np.zeros((n_t, n_phi), dtype=np.float64)

  delta_Q2_total = Q2_edges[-1] - Q2_edges[0]
  delta_nu_total = nu_edges[-1] - nu_edges[0]

  for it in range(n_t):
    delta_t = t_edges[it+1] - t_edges[it]

    for iphi in range(n_phi):
      # Sum over Q2 and nu
      total_var = 0.0
      for iq in range(n_Q2):
        for inu in range(n_nu):
          total_var += var4D_array[it, iq, inu, iphi]

      norm = luminosity * delta_t * delta_Q2_total * delta_nu_total * delta_phi[iphi]
      var_t_phi[it, iphi] = total_var / norm**2

  return var_t_phi

# **********************************
# Get the error for the |t| depedent cross section 
def get_err_t(var_t_phi, delta_phi):
  n_t, n_phi = var_t_phi.shape
  var_t = np.zeros(n_t)

  for it in range(n_t):
    sum = 0.0

    for iphi in range(n_phi):
      sum += var_t_phi[it, iphi] * delta_phi[iphi]**2

    var_t[it] = sum 

  return np.sqrt(var_t)

# **********************************
# Test get_err_S2 function with a small, known-value toy dataset
def toy_test_S_functions():
  # --- Toy array dimensions: [t][Q2][nu][phi] ---
  shape = (2,2,2,2)

  # --- Create toy arrays with small integer counts ---
  real_sum = np.full(shape, 100.0)       # 100 events in each bin
  BH_sum = np.full(shape, 20.0)          # 20 BH events
  lepPi0_sum = np.full(shape, 10.0)      # 10 LEP π0 events
  hepPi0_sum = np.full(shape, 5.0)       # 5 HEP π0 events

  # --- Coefficients ---
  CBH = 1.0
  CPI0_LEP = 1.0
  CPI0_HEP = 1.0
  R_LEPTO = 0.5
  period_idx = 0

  # --- Compute S manually ---
  S_manual = real_sum - CBH*BH_sum - CPI0_LEP*R_LEPTO*lepPi0_sum - CPI0_HEP*(1-R_LEPTO)*hepPi0_sum

  # --- Compute variance manually (Poisson-like) ---
  varS_manual = real_sum + (CBH**2)*BH_sum + (CPI0_LEP*R_LEPTO)**2*lepPi0_sum + (CPI0_HEP*(1-R_LEPTO))**2*hepPi0_sum

  # --- Now use get_S and get_err_S2 (simulating pipeline) ---
  # Wrap arrays in a dict to mimic what get_err_S2 expects
  D_array = real_sum
  B_array = BH_sum
  L_array = lepPi0_sum
  H_array = hepPi0_sum

  # Use the same function your pipeline uses (same mechanics as original)
  def get_err_S2_toy(D_array, B_array, L_array, H_array, period_idx):
    t1 = D_array
    t2 = (CBH)**2 * B_array
    t3 = (CPI0_LEP * R_LEPTO)**2 * L_array
    t4 = (CPI0_HEP * (1-R_LEPTO))**2 * H_array
    return t1 + t2 + t3 + t4

  # Compute variance using the "pipeline" function
  varS_pipeline = get_err_S2_toy(D_array, B_array, L_array, H_array, period_idx)

  # --- Compare ---
  print("=== Toy test of S-term variance ===")
  print("Manual variance:\n", varS_manual)
  print("Pipeline variance:\n", varS_pipeline)
  print("Difference:\n", varS_pipeline - varS_manual)
  print("Relative difference:\n", (varS_pipeline - varS_manual)/varS_manual)

  # --- Check a single bin relative error ---
  idx = (0,0,0,0)
  print("\nExample single bin comparison:")
  print("Manual sqrt(varS):", np.sqrt(varS_manual[idx]))
  print("Pipeline sqrt(varS):", np.sqrt(varS_pipeline[idx]))
  print("Manual relative error:", np.sqrt(varS_manual[idx])/S_manual[idx])
  print("Pipeline relative error:", np.sqrt(varS_pipeline[idx])/S_manual[idx])


# *******************************************************************
# *                   *** T-SLOPE TEST ***                          *
# *******************************************************************
# Get the slope for each period (averaged over both beam charges)
# This value is a used for a spot check - does not follow the true procedure 
def dvcs_slope_period_test(t_edges, sigma_muPlus, sigma_muMinus, err_muPlus, err_muMinus, period):
  # Compute bin centers
  t_edges = np.array(t_edges, dtype=np.float64)
  t_centers = 0.5 * (t_edges[1:] + t_edges[:-1])
  
  # Weighted average over μ+ and μ-
  sigma_avg = (sigma_muPlus / err_muPlus**2 + sigma_muMinus / err_muMinus**2) / \
              (1 / err_muPlus**2 + 1 / err_muMinus**2)
  err_avg = 1 / np.sqrt(1 / err_muPlus**2 + 1 / err_muMinus**2)
  print(sigma_avg)
  print(err_avg)  

  # Remove bad bins (important for stability)
  mask = (err_avg > 0) & np.isfinite(sigma_avg)
  t_fit = t_centers[mask]
  sigma_fit = sigma_avg[mask]
  err_fit = err_avg[mask]
  
  # Model
  def model(t, sigma0, B):
      return sigma0 * np.exp(-B * t)
  
  # Initial guess
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
  plt.errorbar(t_fit, sigma_fit, yerr=err_fit, fmt='o')

  # Fit curve
  t_smooth = np.linspace(min(t_fit), max(t_fit), 200)
  plt.plot(t_smooth, model(t_smooth, sigma0_fit, B_fit), color='black')

  plt.xlabel("|t| (GeV$^2$)")
  plt.ylabel("dsigma/dt (nb/GeV$^2$)")
  plt.title("DVCS t-slope fit")

  # Fit result text
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
  plt.yscale('log')
  plt.savefig(f"dvcs_t_fit{period}.png", dpi=300)
  plt.close()

  #print(f"B = {B_fit:.3f} ± {B_err:.3f} GeV^-2")
  print(f"Saved plot: dvcs_t_fit_{period}.png")
  
  return B_fit, B_err, sigma0_fit, sigma0_err


# *******************************************************************
# *                     *** MAIN ***                                *
# *******************************************************************
# **********************************
# Main function - exclude or include functions here 
def main(): 
  # Testing with toy sample
  #toy_test_S_functions()

  # Dictionary to store the results 
  dvcs_results = {}

  # Loop over the periods 
  debug_main = False
  for idx, period in enumerate(PERIODS[:6]):
    print(idx, period)
    # ***********************************************
    # *          *** 4D ACCEPTANCE ***              * 
    # ***********************************************
    
    gen_muPlus, gen_muMinus, gen_muPlus_sq, gen_muMinus_sq = fill_weights(data_type="gen", period=period)
    rec_muPlus, rec_muMinus, rec_muPlus_sq, rec_muMinus_sq = fill_weights(data_type="rec", period=period)
    acc_muPlus, acc_muMinus = compute_acceptance(rec_muPlus, rec_muMinus, gen_muPlus, gen_muMinus)

    if debug_main:
      test_acceptance(rec_muPlus, rec_muMinus, rec_muPlus_sq, rec_muMinus_sq, gen_muPlus, 
                      gen_muMinus, gen_muPlus_sq, gen_muMinus_sq)

    # ***********************************************
    # *         *** 4D CROSS SECTION ***            * 
    # ***********************************************
    """ """
    # Get the real data sum term in the cross section 
    real_ijkl_muPlus  = np.zeros((4,4,4,8), dtype=np.float64)
    real_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    unweighted_sum(real_ijkl_muPlus, real_ijkl_muMinus, real_files, data="real", period=period)
    if debug_main: 
      print(f"Total real μ+ events: {np.sum(real_ijkl_muPlus)}")
      print(f"Total real μ- events: {np.sum(real_ijkl_muMinus)}")
    
    # Get the BH sum term in the cross section 
    BH_ijkl_muPlus  = np.zeros((4,4,4,8), dtype=np.float64)
    BH_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    weighted_sum(BH_ijkl_muPlus, BH_ijkl_muMinus, hepBH_files, data="hepBH", period=period)
    if debug_main:
      print(f"Total BH μ+ events: {np.sum(BH_ijkl_muPlus)}")
      print(f"Total BH μ- events: {np.sum(BH_ijkl_muMinus)}")
     
    # Get the inclusive inv. pi0 sum term in the cross section
    lepPi0_ijkl_muPlus  = np.zeros((4,4,4,8), dtype=np.float64)
    lepPi0_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    unweighted_sum(lepPi0_ijkl_muPlus, lepPi0_ijkl_muMinus, lepPi0_files, data="lepPi0", period=period)
    if debug_main:
      print(f"Total inc. pi0 μ+ events: {np.sum(lepPi0_ijkl_muPlus)}")
      print(f"Total inc. pi0 μ- events: {np.sum(lepPi0_ijkl_muMinus)}")

    # Get the eexclusive inv. pi0 sum term in the cross section 
    hepPi0_ijkl_muPlus  = np.zeros((4,4,4,8), dtype=np.float64)
    hepPi0_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    weighted_sum(hepPi0_ijkl_muPlus, hepPi0_ijkl_muMinus, hepPi0_files, data="hepPi0", period=period)
    if debug_main:
      print(f"Total exc. pi0 μ+ events: {np.sum(hepPi0_ijkl_muPlus)}")
      print(f"Total exc. pi0 μ- events: {np.sum(hepPi0_ijkl_muMinus)}")

    # Acceptance corrected counts
    sum_ijkl_muPlus = get_S(real_ijkl_muPlus, BH_ijkl_muPlus, lepPi0_ijkl_muPlus, 
                                  hepPi0_ijkl_muPlus, period_idx=idx, charge="muPlus") 
    Ncorr_ijkl_muPlus = (1/acc_muPlus) * sum_ijkl_muPlus

    sum_ijkl_muMinus = get_S(real_ijkl_muMinus, BH_ijkl_muMinus, lepPi0_ijkl_muMinus, 
                                  hepPi0_ijkl_muMinus, period_idx=idx, charge="muMinus") 
    Ncorr_ijkl_muMinus = (1/acc_muMinus) * sum_ijkl_muMinus

    # ***********************************************
    # *     *** T-DEPDENDENT CROSS SECTION ***      * 
    # ***********************************************
    # mu+ t-dependent cross section 
    sigma_t_phi_muPlus = compute_sigma_t_phi(Ncorr_ijkl_muPlus, LUMINOSITY_MUPLUS[idx], 
                                             t_edges, delta_phi, Q2_edges, nu_edges)
    sigma_t_muPlus = integrate_over_phi(sigma_t_phi_muPlus, delta_phi)
    sigma_t_muPlus /= 1e-33 
    print(period, "mu+ dsigma/dt (nb/GeV²):", sigma_t_muPlus)

    # mu- t-dependent cross section 
    sigma_t_phi_muMinus = compute_sigma_t_phi(Ncorr_ijkl_muMinus, LUMINOSITY_MUMINUS[idx], 
                                              t_edges, delta_phi, Q2_edges, nu_edges)
    sigma_t_muMinus = integrate_over_phi(sigma_t_phi_muMinus, delta_phi)
    sigma_t_muMinus /= 1e-33 
    print(period, "mu- dsigma/dt (nb/GeV²):", sigma_t_muMinus)
    
    # ************************************************
    # *       *** 4D CROSS SECTION ERROR ***         * 
    # ************************************************
    # Get the acceptance error 
    varAcc_muPlus, varAcc_muMinus, errAcc_muPlus, errAcc_muMinus = compute_acceptance_error(rec_muPlus_sq, rec_muMinus_sq, gen_muPlus_sq, gen_muMinus_sq,
                                                      rec_muPlus, rec_muMinus, gen_muPlus, gen_muMinus) 
    
    # Get the sum term error 
    dD_ijkl_muPlus = np.zeros((4,4,4,8), dtype=np.float64)
    dD_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    unweighted_err_sum(dD_ijkl_muPlus, dD_ijkl_muMinus, real_files, data="real", period=period)
    if debug_main: 
      print("mu+ varD min/max:", np.min(dD_ijkl_muPlus), np.max(dD_ijkl_muPlus))
      print("mu- varD min/max:", np.min(dD_ijkl_muMinus), np.max(dD_ijkl_muMinus))
    
    dB_ijkl_muPlus = np.zeros((4,4,4,8), dtype=np.float64)
    dB_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    weighted_err_sum(dB_ijkl_muPlus, dB_ijkl_muMinus, hepBH_files, data="hepBH", period=period)
    if debug_main:
      print("mu+ varB min/max:", np.min(dB_ijkl_muPlus), np.max(dB_ijkl_muPlus))
      print("mu- varB min/max:", np.min(dB_ijkl_muMinus), np.max(dB_ijkl_muMinus))

    dL_ijkl_muPlus = np.zeros((4,4,4,8), dtype=np.float64)
    dL_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    unweighted_err_sum(dL_ijkl_muPlus, dL_ijkl_muMinus, lepPi0_files, data="lepPi0", period=period)
    if debug_main:
      print("mu+ varL min/max:", np.min(dL_ijkl_muPlus), np.max(dL_ijkl_muPlus))
      print("mu- varL min/max:", np.min(dL_ijkl_muMinus), np.max(dL_ijkl_muMinus))

    dH_ijkl_muPlus = np.zeros((4,4,4,8), dtype=np.float64)
    dH_ijkl_muMinus = np.zeros((4,4,4,8), dtype=np.float64)
    weighted_err_sum(dH_ijkl_muPlus, dH_ijkl_muMinus, hepPi0_files, data="hepPi0", period=period)
    if debug_main:
      print("mu+ varH min/max:", np.min(dH_ijkl_muPlus), np.max(dH_ijkl_muPlus))
      print("mu- varH min/max:", np.min(dH_ijkl_muMinus), np.max(dH_ijkl_muMinus))
    
    dS_ijkl_muPlus = get_err_S2(dD_ijkl_muPlus, dB_ijkl_muPlus, dL_ijkl_muPlus, dH_ijkl_muPlus, period_idx=idx, charge="muPlus")
    dS_ijkl_muMinus = get_err_S2(dD_ijkl_muMinus, dB_ijkl_muMinus, dL_ijkl_muMinus, dH_ijkl_muMinus, period_idx=idx, charge="muMinus")
    if debug_main:
      print("mu+ varS min/max:", np.min(dS_ijkl_muPlus), np.max(dS_ijkl_muPlus))
      print("mu- varS min/max:", np.min(dS_ijkl_muMinus), np.max(dS_ijkl_muMinus))

    if debug_main: # detailed error calculation check 
      print("Error diagnostics (muPlus only):")
      print("S_ijkl: min =", np.min(sum_ijkl_muPlus), "max =", np.max(sum_ijkl_muPlus))
      print("A_ijkl: min =", np.min(acc_muPlus), "max =", np.max(acc_muPlus))
      print("varS_ijkl: min =", np.min(dS_ijkl_muPlus), "max =", np.max(dS_ijkl_muPlus))
      print("varA_ijkl: min =", np.min(varAcc_muPlus), "max =", np.max(varAcc_muPlus))
      print("S/A ratio: min =", np.min(sum_ijkl_muPlus/acc_muPlus), "max =", np.max(sum_ijkl_muPlus/acc_muPlus))
      print("varS/S^2: min =", np.min(dS_ijkl_muPlus/(sum_ijkl_muPlus**2)), "max =", np.max(dS_ijkl_muPlus/(sum_ijkl_muPlus**2)))
      print("varA/A^2: min =", np.min(varAcc_muPlus/(acc_muPlus**2)), "max =", np.max(varAcc_muPlus/(acc_muPlus**2)))
      print("(S/A)^2 ratio: min =", np.min((sum_ijkl_muPlus**2)/(acc_muPlus**2)), "max =", np.max((sum_ijkl_muPlus**2)/(acc_muPlus**2)))

    # ************************************************
    # *   *** T-DEPDENDENT CROSS SECTION ERROR ***   * 
    # ************************************************
    # mu+ t-dependent cross section  error 
    err_ijkl_muPlus = compute_4D_error(sum_ijkl_muPlus, acc_muPlus, dS_ijkl_muPlus, varAcc_muPlus)
    err_t_phi_muPlus = get_err_t_phi(err_ijkl_muPlus, LUMINOSITY_MUPLUS[idx], t_edges, delta_phi, Q2_edges, nu_edges)
    err_t_muPlus = get_err_t(err_t_phi_muPlus, delta_phi)
    err_t_muPlus /= 1e-33
    print(period, "mu+ error (nb/GeV²):", err_t_muPlus)

    # mu- t-dependent cross section error
    err_ijkl_muMinus = compute_4D_error(sum_ijkl_muMinus, acc_muMinus, dS_ijkl_muMinus, varAcc_muMinus)
    err_t_phi_muMinus = get_err_t_phi(err_ijkl_muMinus, LUMINOSITY_MUMINUS[idx], t_edges, delta_phi, Q2_edges, nu_edges)
    err_t_muMinus = get_err_t(err_t_phi_muMinus, delta_phi)
    err_t_muMinus /= 1e-33
    print(period, "mu- error (nb/GeV²):", err_t_muMinus)
    """ """

    # Spot check to see if slopes seem reasonable - not actual fitting method but good approx. for a check
    # If the numbers seem unreasonable there is a problem somewhere 
    if debug_main:
      B, B_err, sigma0, sigma0_err = dvcs_slope_period_test(t_edges, sigma_t_muPlus, sigma_t_muMinus, err_t_muPlus, err_t_muMinus, period)
      print(f"{period} Slope Test: B = {B:.3f} ± {B_err:.3f} GeV^-2, sigma0 = {sigma0:.2f} ± {sigma0_err:.2f} nb/GeV^2")

    # Save results to the dictionary 
    dvcs_results[period] = {
      "sigma_muPlus": sigma_t_muPlus.copy(),
      "sigma_muMinus": sigma_t_muMinus.copy(),
      "err_muPlus": err_t_muPlus.copy(),
      "err_muMinus": err_t_muMinus.copy()
    }
    """ """
  # Save dictionary to output file for later use 
  with open("dvcs_xSection_results.pkl", "wb") as f:
    pickle.dump(dvcs_results, f)
  print("Results written to 'dvcs_xSection_results.pkl'")

if __name__ == "__main__":
  main()