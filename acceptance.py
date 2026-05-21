import os
import ROOT
from ROOT import TLorentzVector
import math
import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal, getcontext
from matplotlib import gridspec
import datetime  
import pickle

# TODO: Use all data instead of just P04  - just one period is used right now for testing 

# *******************************************************************
# *                    *** DATA PREP ***                            *
# *******************************************************************
# Generated MC data 
gen_dir = "/eos/user/g/gkainth/BH/"
gen_files = [os.path.join(gen_dir, "gen_P04_muPlus.root"), os.path.join(gen_dir, "gen_P04_muMinus.root")]

""" gen_files = [os.path.join(gen_dir, "gen_P04_muPlus.root"), os.path.join(gen_dir, "gen_P04_muMinus.root"),
             os.path.join(gen_dir, "gen_P05_muPlus.root"), os.path.join(gen_dir, "gen_P05_muMinus.root"),
             os.path.join(gen_dir, "gen_P06_muPlus.root"), os.path.join(gen_dir, "gen_P06_muMinus.root"),
             os.path.join(gen_dir, "gen_P07_muPlus.root"), os.path.join(gen_dir, "gen_P07_muMinus.root"),
             os.path.join(gen_dir, "gen_P08_muPlus.root"), os.path.join(gen_dir, "gen_P08_muMinus.root"),
             os.path.join(gen_dir, "gen_P09_muPlus.root"), os.path.join(gen_dir, "gen_P09_muMinus.root")] """

tree_gen_muPlus = ROOT.TChain("USR970_gen")
tree_gen_muMinus = ROOT.TChain("USR970_gen")

for f in gen_files:
  if "muPlus" in f:
    tree_gen_muPlus.Add(f)
  elif "muMinus" in f:
    tree_gen_muMinus.Add(f)

print("Total Entries Generated: mu+:", tree_gen_muPlus.GetEntries(), ", mu-:", tree_gen_muMinus.GetEntries())

# **********************************
# Reconstructed MC data after full DVCS selection
rec_dir = "/eos/user/g/gkainth/BH/"
rec_files = [os.path.join(rec_dir, "filtered_P04_muPlus.root"), os.path.join(rec_dir, "filtered_P04_muMinus.root")]

""" rec_files = [os.path.join(rec_dir, "filtered_P04_muPlus.root"), os.path.join(rec_dir, "filtered_P04_muMinus.root"),
             os.path.join(rec_dir, "filtered_P05_muPlus.root"), os.path.join(rec_dir, "filtered_P05_muMinus.root"),
             os.path.join(rec_dir, "filtered_P06_muPlus.root"), os.path.join(rec_dir, "filtered_P06_muMinus.root"),
             os.path.join(rec_dir, "filtered_P07_muPlus.root"), os.path.join(rec_dir, "filtered_P07_muMinus.root"),
             os.path.join(rec_dir, "filtered_P08_muPlus.root"), os.path.join(rec_dir, "filtered_P08_muMinus.root"),
             os.path.join(rec_dir, "filtered_P09_muPlus.root"), os.path.join(rec_dir, "filtered_P09_muMinus.root")] """

tree_rec_muPlus = ROOT.TChain("USR970_filtered")
tree_rec_muMinus = ROOT.TChain("USR970_filtered")

for f in rec_files:
  if "muPlus" in f:
    tree_rec_muPlus.Add(f)
  elif "muMinus" in f:
    tree_rec_muMinus.Add(f)

print("Total Entries Reconstructed: mu+:", tree_rec_muPlus.GetEntries(), ", mu-:", tree_rec_muMinus.GetEntries())


# *******************************************************************
# *              *** BINNED PHASE SPACE PLOTS ***                   *
# *******************************************************************
# Function for drawing the bin lines 
def draw_bin_lines_2D(hist, x_edges, y_edges):
  lines = []  # store lines to keep them alive

  # Vertical lines (x-axis bin edges)
  for x in x_edges:
    line = ROOT.TLine(x, y_edges[0], x, y_edges[-1])
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(1)
    line.Draw("same")
    lines.append(line)

  # Horizontal lines (y-axis bin edges)
  for y in y_edges:
    line = ROOT.TLine(x_edges[0], y, x_edges[-1], y)
    line.SetLineColor(ROOT.kBlack)
    line.SetLineWidth(1)
    line.Draw("same")
    lines.append(line)

  return lines

# **********************************
# Function for making default plots 
def default_plots():
  # Define histograms
  hMC_Q2nu = ROOT.TH2F("hMC_Q2nu", "Q^{2}-#nu Distribution - Reconstructed MC; Q^{2} [(GeV/c)^{2}]; #nu [GeV]", 100, 0, 11, 100, 0, 35)
  hMC_Q2t = ROOT.TH2F("hMC_Q2t", "Q^{2}-|t| Distribution - Reconstructed MC; Q^{2} [(GeV/c)^{2}]; |t| [(GeV/c)^{2}]", 100, 0, 11, 100, 0, 1)

  # Loop over events and fill histograms
  for event in tree_rec_muPlus:
    Q2 = event.Q2_fit
    nu = event.nu_fit
    t = event.t_fit

    hMC_Q2nu.Fill(Q2, nu)
    hMC_Q2t.Fill(Q2, abs(t))

  for event in tree_rec_muMinus:
    Q2 = event.Q2_fit
    nu = event.nu_fit
    t = event.t_fit

    hMC_Q2nu.Fill(Q2, nu)
    hMC_Q2t.Fill(Q2, abs(t))

  # Create canvases to draw and save the histograms
  c1 = ROOT.TCanvas("c1", "Q2 vs nu", 800, 600)
  hMC_Q2nu.Draw("COLZ")
  lines1 = draw_bin_lines_2D(hMC_Q2nu, Q2_edges, nu_edges)

  c1.Update()  # make sure the stats box is created
  st1 = hMC_Q2nu.GetListOfFunctions().FindObject("stats")
  if st1:
    st1.SetX1NDC(0.72)  # left edge (moved a bit left for more width)
    st1.SetX2NDC(0.93)  # right edge
    st1.SetY1NDC(0.70)  # bottom edge (a little lower for extra height)
    st1.SetY2NDC(0.89)  # top edge
    st1.SetTextSize(0.03)
    st1.SetFillColorAlpha(ROOT.kWhite, 0.6)
    st1.SetBorderSize(1)

  hMC_Q2nu.GetZaxis().SetLabelSize(0.025)
  hMC_Q2nu.GetZaxis().SetTitleOffset(1.2)
  c1.SetRightMargin(0.15)
  c1.SaveAs("Q2_vs_nu.png")

  c2 = ROOT.TCanvas("c2", "Q2 vs |t|", 800, 600)
  hMC_Q2t.Draw("COLZ")
  lines2 = draw_bin_lines_2D(hMC_Q2t, Q2_edges, t_edges)

  c2.Update()  # make sure the stats box is created
  st2 = hMC_Q2t.GetListOfFunctions().FindObject("stats")
  if st2:
    st2.SetX1NDC(0.72)  # left edge (moved a bit left for more width)
    st2.SetX2NDC(0.93)  # right edge
    st2.SetY1NDC(0.70)  # bottom edge (a little lower for extra height)
    st2.SetY2NDC(0.89)  # top edge
    st2.SetTextSize(0.03)
    st2.SetFillColorAlpha(ROOT.kWhite, 0.5)
    st2.SetBorderSize(1)

  hMC_Q2t.GetZaxis().SetLabelSize(0.025)
  hMC_Q2t.GetZaxis().SetTitleOffset(1.2)
  c2.SetRightMargin(0.15)
  c2.SaveAs("Q2_vs_t.png")


# *******************************************************************
# *              *** 3D ACCEPTANCE FUNCTIONS ***                    *
# *******************************************************************
# Define the bin edges/bins for wide acceptance calculation 
# nu: 11 bins of width 2 GeV between 10 and 32 GeV
nu_edges = np.linspace(10, 32, 12)
#nu_edges = np.linspace(10, 32, 5)
n_nu_bins = len(nu_edges) - 1
# Q2: 9 bins of width 1 (GeV/c)^2 between 1 and 10
Q2_edges = np.linspace(1, 10, 10)
n_Q2_bins = len(Q2_edges) - 1
# |t|: 4 bins -> each bin should have roughly the same no. of events 
t_edges = [0.08, 0.136, 0.219, 0.36, 0.64]
n_t_bins = len(t_edges) - 1
# phi: 8 bins of width pi/4 rad between -pi and pi
phi_edges = np.linspace(-np.pi, np.pi, 9)
n_phi_bins = len(phi_edges) - 1

# **********************************
# Function for initializing 3D weight arrays for acceptance and error 
# Returns a dictionary of arrays 
def init_acceptance_arrays(shape):
  return {
    "rec_plus":  np.zeros(shape, dtype=np.float64),
    "rec_minus": np.zeros(shape, dtype=np.float64),
    "gen_plus":  np.zeros(shape, dtype=np.float64),
    "gen_minus": np.zeros(shape, dtype=np.float64),

    "rec2_plus":  np.zeros(shape, dtype=np.float64),
    "rec2_minus": np.zeros(shape, dtype=np.float64),
    "gen2_plus":  np.zeros(shape, dtype=np.float64),
    "gen2_minus": np.zeros(shape, dtype=np.float64)
  }

# **********************************
# Function for filling the 3D weight arrays
def fill_weights_3d(tree, charge, data_type, arrays, project="t", Q2_edges=Q2_edges, 
                 nu_edges=nu_edges, phi_edges=phi_edges, t_edges=t_edges):
  for event in tree:
    weight = np.float64(event.weight_DVCS)

    # Extract variables
    if data_type == "gen":
      Q2  = np.float64(event.Q2_gen)
      nu  = np.float64(event.nu_gen)
      t   = np.float64(event.t_gen)
      phi = np.float64(event.phi_gg_gen)

    elif data_type == "rec":
      Q2  = np.float64(event.Q2_fit)
      nu  = np.float64(event.nu_fit)
      t   = np.float64(event.t_fit)
      phi = np.float64(event.phi_gg_fit)

    else:
      raise ValueError("data_type must be 'gen' or 'rec'")

    # Q2 bin
    Q2_bin = np.searchsorted(Q2_edges, Q2) - 1
    if Q2_bin < 0 or Q2_bin >= len(Q2_edges) - 1:
      continue

    # nu bin
    nu_bin = np.searchsorted(nu_edges, nu) - 1
    if nu_bin < 0 or nu_bin >= len(nu_edges) - 1:
      continue

    if project == "t":
      # Bin phi
      phi = ((phi + np.pi) % (2*np.pi)) - np.pi
      phi_bin = np.searchsorted(phi_edges, phi) - 1
      if phi_bin < 0 or phi_bin >= len(phi_edges) - 1:
        continue
      #print("Q2:", Q2_bin, Q2)
      #print("nu:", nu_bin, nu)
      #print("phi:", phi_bin, phi)
      # Axis order: [Q2][nu][phi]
      i1, i2, i3 = Q2_bin, nu_bin, phi_bin
      
    elif project == "phi":
      # Bin |t|
      t_abs = abs(t)
      t_bin = np.searchsorted(t_edges, t_abs) - 1
      if t_bin < 0 or t_bin >= len(t_edges) - 1:
        continue
      #print("Q2:", Q2_bin, Q2)
      #print("nu:", nu_bin, nu)
      #print("t:", t_bin, t_abs)
      # Axis order: [t][Q2][nu]
      i1, i2, i3 = t_bin, Q2_bin, nu_bin

    else:
      raise ValueError("project must be 't' or 'phi'")

    # Fill the appropriate arrays
    if data_type == "rec":
      if charge == "muPlus":
        arrays["rec_plus"][i1, i2, i3] += weight
        arrays["rec2_plus"][i1, i2, i3] += weight ** 2
      elif charge == "muMinus":
        arrays["rec_minus"][i1, i2, i3] += weight
        arrays["rec2_minus"][i1, i2, i3] += weight ** 2
    elif data_type == "gen":
      if charge == "muPlus":
        arrays["gen_plus"][i1, i2, i3] += weight
        arrays["gen2_plus"][i1, i2, i3] += weight ** 2
      elif charge == "muMinus":
        arrays["gen_minus"][i1, i2, i3] += weight
        arrays["gen2_minus"][i1, i2, i3] += weight ** 2

# **********************************
# Function for calculating the error (returns result for mu+ and mu-)
def acc_error_3d(arrays, project="t", Q2_edges=Q2_edges,
              nu_edges=nu_edges, phi_edges=phi_edges, t_edges=t_edges):

  # Select shape based on projection
  if project == "t":
    shape = (len(Q2_edges) - 1,
             len(nu_edges) - 1,
             len(phi_edges) - 1)

  elif project == "phi":
    shape = (len(t_edges) - 1,
             len(Q2_edges) - 1,
             len(nu_edges) - 1)

  else:
    raise ValueError("project must be 't' or 'phi'")

  # Output arrays
  acc_err_muPlus = np.zeros(shape, dtype=np.float64)
  acc_err_muMinus = np.zeros(shape, dtype=np.float64)

  # Loop over bins
  for idx in np.ndindex(shape):
    # mu+
    r = arrays["rec_plus"][idx]
    g = arrays["gen_plus"][idx]
    r2 = arrays["rec2_plus"][idx]
    g2 = arrays["gen2_plus"][idx]

    if g == 0:
      acc_err_muPlus[idx] = 0.0
    else:
      term1 = (1.0 / g) ** 2 * r2
      term2 = (r / (g ** 2)) ** 2 * g2
      acc_err_muPlus[idx] = np.sqrt(term1 + term2)

    # mu-
    r = arrays["rec_minus"][idx]
    g = arrays["gen_minus"][idx]
    r2 = arrays["rec2_minus"][idx]
    g2 = arrays["gen2_minus"][idx]

    if g == 0:
      acc_err_muMinus[idx] = 0.0
    else:
      term1 = (1.0 / g) ** 2 * r2
      term2 = (r / (g ** 2)) ** 2 * g2
      acc_err_muMinus[idx] = np.sqrt(term1 + term2)

  return acc_err_muPlus, acc_err_muMinus

# **********************************
# Plot the dsitrubtion of the events by bin (projected in t)
def plot_acceptance_project_t(acc_muPlus, acc_muMinus, acc_err_muPlus, acc_err_muMinus):
  phi_bin_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])
  n_Q2_bins, n_nu_bins, n_phi_bins = acc_muPlus.shape

  fig, axes = plt.subplots(nrows=n_nu_bins, ncols=n_Q2_bins, figsize=(18, 22), sharex=True, sharey=True)

  for i in range(n_nu_bins):    # nu bins (y-axis)
    for j in range(n_Q2_bins):  # Q2 bins (x-axis)
      ax = axes[i, j]

      y_muPlus = acc_muPlus[j, i]
      yerr_muPlus = acc_err_muPlus[j, i]
      y_muMinus = acc_muMinus[j, i]
      yerr_muMinus = acc_err_muMinus[j, i]

      ax.errorbar(phi_bin_centers, y_muPlus, yerr=yerr_muPlus, fmt='o', markerfacecolor='none', 
                  markeredgecolor='red', markersize=5, ecolor='red', label='μ⁺' if i == 0 and j == 0 else "")
      ax.errorbar(phi_bin_centers, y_muMinus, yerr=yerr_muMinus, fmt='o', markerfacecolor='none', 
                  markeredgecolor='black', markersize=5, ecolor='black', label='μ⁻' if i == 0 and j == 0 else "")
      
      ax.set_ylim(0, 0.8)
      ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
      ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"], fontsize=10)
      ax.grid(True, linestyle='--', linewidth=0.5)
      ax.axhline(0, color='gray', linewidth=0.5)

  # Legend
  handles, labels = axes[0, 0].get_legend_handles_labels()
  fig.legend(handles, labels, loc='upper right', fontsize=14, markerscale=1.5)

  # Create a new set of axes for the phi and acceptance scale at the top right corner
  phi_axis = fig.add_axes([0.817, 0.92, 0.076, 0.034])  # [left, bottom, width, height]
  acc_axis = fig.add_axes([0.8, 0.877, 0.11, 0.064])

  # Setup phi axis
  phi_axis.set_xlim(-np.pi, np.pi)
  phi_axis.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
  phi_axis.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"], fontsize=12)
  phi_axis.set_yticks([])
  phi_axis.yaxis.set_visible(False)
  phi_axis.tick_params(axis='x', direction='in', length=5, top=True, bottom=False)
  phi_axis.xaxis.set_label_position('top')
  phi_axis.set_xlabel(r"$\phi_{\gamma^*\gamma}$ [rad]", fontsize=14, labelpad=20)
  phi_axis.xaxis.tick_top()
  phi_axis.patch.set_facecolor('none')
  for name, spine in phi_axis.spines.items():
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

  # Global axes for Q2 and nu 
  Q2_axis = fig.add_axes([0.116, 0.09, 0.778, 0.035])
  nu_axis = fig.add_axes([0.085, 0.11, 0.035, 0.835])

  # Q2 axis setup
  Q2_axis.set_xlim(Q2_edges[0], Q2_edges[-1])
  Q2_axis.set_xticks(Q2_edges)
  Q2_axis.set_xticklabels([f"{int(a)}" for a in Q2_edges], fontsize=12)
  Q2_axis.set_yticks([])
  Q2_axis.xaxis.tick_bottom()
  Q2_axis.set_xlabel(r"$Q^2$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
  Q2_axis.patch.set_facecolor('none')
  for name, spine in Q2_axis.spines.items():
    spine.set_visible(name == 'bottom')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  # nu axis setup
  nu_axis.set_ylim(nu_edges[-1], nu_edges[0])
  nu_axis.set_yticks(nu_edges)
  nu_axis.set_yticklabels([f"{int(a)}" for a in nu_edges], fontsize=12)
  nu_axis.set_xticks([])
  nu_axis.yaxis.tick_left()
  nu_axis.set_ylabel(r"$\nu$ [GeV]", fontsize=14, labelpad=20)
  nu_axis.patch.set_facecolor('none')
  for name, spine in nu_axis.spines.items():
    spine.set_visible(name == 'left')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  plt.tight_layout(rect=[0.10, 0.10, 0.9, 0.95])
  plt.savefig("acceptance_projected_t.png", dpi=300)


# **********************************
# Plot the dsitrubtion of the events by bin (projected in phi)
def plot_acceptance_project_phi(acc_muPlus, acc_muMinus, acc_err_muPlus, acc_err_muMinus):
  nu_bin_centers = 0.5 * (nu_edges[:-1] + nu_edges[1:])
  n_t_bins, n_Q2_bins, n_nu_bins = acc_muPlus.shape

  fig, axes = plt.subplots(nrows=n_Q2_bins, ncols=n_t_bins, figsize=(18, 22), sharex=True, sharey=True)

  for i in range(n_Q2_bins):  # rows (Q2 -> y values)    
    for j in range(n_t_bins):  # columns (t -> x)   
      ax = axes[i, j]

      # acceptance_muPlus and acceptance_muMinus shape: (t, Q2, nu)
      y_muPlus = acc_muPlus[j, i]
      yerr_muPlus = acc_err_muPlus[j, i]
      y_muMinus = acc_muMinus[j, i]
      yerr_muMinus = acc_err_muMinus[j, i]

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
  plt.savefig("acceptance_projected_phi.png", dpi=300)


# *******************************************************************
# *              *** 4D ACCEPTANCE FUNCTIONS ***                    *
# *******************************************************************
# Tighten the binning in Q2 and nu 
# Q2 and nu (tight bins)
Q2_edges_tight = np.linspace(1, 5, 5)   # 4 bins
nu_edges_tight = np.linspace(10, 32, 5) # 4 bins

# **********************************
# Fill the weights in the 4D arrays 
def fill_weights_4d(tree, charge, data_type, arrays, Q2_edges=Q2_edges_tight, nu_edges=nu_edges_tight, 
                    t_edges=t_edges, phi_edges=phi_edges, debug=False):

  # Debugging: track how many events pass each stage
  total_events = 0
  filled_events = 0

  for event in tree:
    total_events += 1

    weight = np.float64(event.weight_DVCS)

    if data_type == "gen":
      Q2  = np.float64(event.Q2_gen)
      nu  = np.float64(event.nu_gen)
      t   = np.float64(event.t_gen)
      phi = np.float64(event.phi_gg_gen)

    elif data_type == "rec":
      Q2  = np.float64(event.Q2_fit)
      nu  = np.float64(event.nu_fit)
      t   = np.float64(event.t_fit)
      phi = np.float64(event.phi_gg_fit)

    else:
      raise ValueError("data_type must be 'gen' or 'rec'")

    # Bin Q2
    Q2_bin = np.searchsorted(Q2_edges, Q2) - 1
    if Q2_bin < 0 or Q2_bin >= len(Q2_edges) - 1:
      continue

    # Bin nu
    nu_bin = np.searchsorted(nu_edges, nu) - 1
    if nu_bin < 0 or nu_bin >= len(nu_edges) - 1:
      continue

    # Bin |t|
    t_abs = abs(t)
    t_bin = np.searchsorted(t_edges, t_abs) - 1
    if t_bin < 0 or t_bin >= len(t_edges) - 1:
      continue

    # Bin phi (wrapped)
    phi = ((phi + np.pi) % (2*np.pi)) - np.pi
    phi_bin = np.searchsorted(phi_edges, phi) - 1
    if phi_bin < 0 or phi_bin >= len(phi_edges) - 1:
      continue

    # Final 4D index
    # Axis order: [Q2][nu][t][phi]
    idx = (Q2_bin, nu_bin, t_bin, phi_bin)

    # Fill arrays
    if data_type == "rec":
      if charge == "muPlus":
        arrays["rec_plus"][idx] += weight
        arrays["rec2_plus"][idx] += weight**2
      elif charge == "muMinus":
        arrays["rec_minus"][idx] += weight
        arrays["rec2_minus"][idx] += weight**2

    elif data_type == "gen":
      if charge == "muPlus":
        arrays["gen_plus"][idx] += weight
        arrays["gen2_plus"][idx] += weight**2
      elif charge == "muMinus":
        arrays["gen_minus"][idx] += weight
        arrays["gen2_minus"][idx] += weight**2

    filled_events += 1

    # Debug prints (optional)
    if debug and filled_events < 10:
      print("Event:", filled_events)
      print("Q2, nu, t, phi:", Q2, nu, t, phi)
      print("Bins:", idx)
      print("---")
      print(f"[{charge} | {data_type}] Total events: {total_events}")
      print(f"[{charge} | {data_type}] Filled events: {filled_events}")


# **********************************
# Function for calculating the error (returns result for mu+ and mu-)
def acc_error_4d(arrays):
  shape = arrays["rec_plus"].shape
  acc_err_muPlus = np.zeros(shape, dtype=np.float64)
  acc_err_muMinus = np.zeros(shape, dtype=np.float64)

  # mu+
  r  = arrays["rec_plus"]
  g  = arrays["gen_plus"]
  r2 = arrays["rec2_plus"]
  g2 = arrays["gen2_plus"]

  mask = g != 0
  term1 = np.zeros(shape)
  term2 = np.zeros(shape)
  term1[mask] = (1.0 / g[mask])**2 * r2[mask]
  term2[mask] = (r[mask] / (g[mask]**2))**2 * g2[mask]
  acc_err_muPlus[mask] = np.sqrt(term1[mask] + term2[mask])

  # mu-
  r  = arrays["rec_minus"]
  g  = arrays["gen_minus"]
  r2 = arrays["rec2_minus"]
  g2 = arrays["gen2_minus"]

  mask = g != 0
  term1 = np.zeros(shape)
  term2 = np.zeros(shape)
  term1[mask] = (1.0 / g[mask])**2 * r2[mask]
  term2[mask] = (r[mask] / (g[mask]**2))**2 * g2[mask]
  acc_err_muMinus[mask] = np.sqrt(term1[mask] + term2[mask])

  return acc_err_muPlus, acc_err_muMinus

# **********************************
#Plot the dsitrubtion of the events by |t| bin (projected over t)
def plot_acceptance_by_tBin_4D(arrays, acc_err_plus, acc_err_minus, idx):
  phi_bin_centers = 0.5 * (phi_edges[:-1] + phi_edges[1:])

  # slice in t
  rec_plus = arrays["rec_plus"][:, :, idx, :]
  gen_plus = arrays["gen_plus"][:, :, idx, :]
  rec_minus = arrays["rec_minus"][:, :, idx, :]
  gen_minus = arrays["gen_minus"][:, :, idx, :]

  # acceptance
  acc_plus = np.zeros_like(rec_plus)
  acc_minus = np.zeros_like(rec_minus)

  mask_plus = gen_plus != 0
  mask_minus = gen_minus != 0

  acc_plus[mask_plus] = rec_plus[mask_plus] / gen_plus[mask_plus]
  acc_minus[mask_minus] = rec_minus[mask_minus] / gen_minus[mask_minus]

  n_Q2_bins, n_nu_bins, n_phi_bins = rec_plus.shape

  fig, axes = plt.subplots(nrows=n_nu_bins, ncols=n_Q2_bins,
                           figsize=(22, 22), sharex=True, sharey=True)

  for i in range(n_nu_bins):
    for j in range(n_Q2_bins):
      # (0,0) is top left corner of plot, need to start from bottom 
      ax = axes[(n_nu_bins - 1 - i), j]

      y_muPlus = acc_plus[j, i, :]
      y_muMinus = acc_minus[j, i, :]

      # errors (same indexing logic as before)
      yerr_muPlus = acc_err_plus[j, i, idx, :]
      yerr_muMinus = acc_err_minus[j, i, idx, :]

      ax.errorbar(phi_bin_centers, y_muPlus, yerr=yerr_muPlus,
                  fmt='o', markerfacecolor='none',
                  markeredgecolor='red', markersize=7, ecolor='red',
                  label='μ⁺' if i == 0 and j == 0 else "")

      ax.errorbar(phi_bin_centers, y_muMinus, yerr=yerr_muMinus,
                  fmt='o', markerfacecolor='none',
                  markeredgecolor='black', markersize=7, ecolor='black',
                  label='μ⁻' if i == 0 and j == 0 else "")

      ax.set_ylim(0, 0.5)
      ax.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
      ax.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0",
                          r"$\pi/2$", r"$\pi$"], fontsize=10)

      ax.grid(True, linestyle='--', linewidth=0.5)
      ax.axhline(0, color='gray', linewidth=0.5)

  handles, labels = axes[0, 0].get_legend_handles_labels()
  fig.legend(handles, labels, loc='upper right', fontsize=14, markerscale=1.5)

    # Create a new set of axes for the phi and acceptance scale at the top right corner
  phi_axis = fig.add_axes([0.71, 0.91, 0.18, 0.05])  # [left, bottom, width, height]
  acc_axis = fig.add_axes([0.86, 0.745, 0.05, 0.195])

  # Setup phi axis
  phi_axis.set_xlim(-np.pi, np.pi)
  phi_axis.set_xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi])
  phi_axis.set_xticklabels([r"$-\pi$", r"$-\pi/2$", "0", r"$\pi/2$", r"$\pi$"], fontsize=12)
  phi_axis.set_yticks([])
  phi_axis.yaxis.set_visible(False)
  phi_axis.tick_params(axis='x', direction='in', length=5, top=True, bottom=False)
  phi_axis.xaxis.set_label_position('top')
  phi_axis.set_xlabel(r"$\phi_{\gamma^*\gamma}$ [rad]", fontsize=14, labelpad=20)
  phi_axis.xaxis.tick_top()
  phi_axis.patch.set_facecolor('none')
  for name, spine in phi_axis.spines.items():
    spine.set_visible(name == 'top')
    if name == 'top':
      spine.set_linewidth(1.0)
      spine.set_color('black')

  # Setup acceptance axis
  acc_axis.set_ylim(0, 0.5)
  acc_axis.set_yticks([0.0, 0.1, 0.2, 0.3, 0.4, 0.5])
  acc_axis.set_yticklabels(["0", "0.1", "0.2", "0.3", "0.4", "0.5"], fontsize=12)
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

  # Global axes for Q2 and nu 
  Q2_axis = fig.add_axes([0.116, 0.09, 0.778, 0.035])
  nu_axis = fig.add_axes([0.085, 0.11, 0.035, 0.835])

  # Q2 axis setup
  Q2_axis.set_xlim(Q2_edges_tight[0], Q2_edges_tight[-1])
  Q2_axis.set_xticks(Q2_edges_tight)
  Q2_axis.set_xticklabels([f"{int(a)}" for a in Q2_edges_tight], fontsize=12)
  Q2_axis.set_yticks([])
  Q2_axis.xaxis.tick_bottom()
  Q2_axis.set_xlabel(r"$Q^2$ [(GeV/c)$^2$]", fontsize=14, labelpad=20)
  Q2_axis.patch.set_facecolor('none')
  for name, spine in Q2_axis.spines.items():
    spine.set_visible(name == 'bottom')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  # nu axis setup
  nu_axis.set_ylim(nu_edges_tight[0], nu_edges_tight[-1])
  nu_axis.set_yticks(nu_edges_tight)
  nu_axis.set_yticklabels([f"{int(a)}" for a in nu_edges_tight], fontsize=12)
  nu_axis.set_xticks([])
  nu_axis.yaxis.tick_left()
  nu_axis.set_ylabel(r"$\nu$ [GeV]", fontsize=14, labelpad=20)
  nu_axis.patch.set_facecolor('none')
  for name, spine in nu_axis.spines.items():
    spine.set_visible(name == 'left')
    spine.set_linewidth(1.0)
    spine.set_color('black')

  plt.tight_layout(rect=[0.10, 0.10, 0.9, 0.95])
  plt.savefig(f"acceptance_integrated_bin{idx+1}.png", dpi=300)



# *******************************************************************
# *                  *** WRAPPER FUNCTIONS ***                      *
# *******************************************************************
def wrapper_project_t():
  shape_int_t = (n_Q2_bins, n_nu_bins, n_phi_bins)
  arrays_int_t = init_acceptance_arrays(shape_int_t)
  
  fill_weights_3d(tree=tree_gen_muPlus, charge="muPlus", data_type="gen", arrays=arrays_int_t, project="t")
  fill_weights_3d(tree=tree_gen_muMinus, charge="muMinus", data_type="gen", arrays=arrays_int_t, project="t")
  fill_weights_3d(tree=tree_rec_muPlus, charge="muPlus", data_type="rec", arrays=arrays_int_t, project="t")
  fill_weights_3d(tree=tree_rec_muMinus, charge="muMinus", data_type="rec", arrays=arrays_int_t, project="t")
  # acceptance mu+ 
  acc_muPlus = np.zeros_like(arrays_int_t["rec_plus"])
  mask_muPlus = arrays_int_t["gen_plus"] != 0
  acc_muPlus[mask_muPlus] = ( arrays_int_t["rec_plus"][mask_muPlus] / arrays_int_t["gen_plus"][mask_muPlus]) 
  # acceptance mu- 
  acc_muMinus = np.zeros_like(arrays_int_t["rec_minus"])
  mask_muMinus = arrays_int_t["gen_minus"] != 0
  acc_muMinus[mask_muMinus] = ( arrays_int_t["rec_minus"][mask_muMinus] / arrays_int_t["gen_minus"][mask_muMinus]) 
  print("mu+ acceptance min:", np.min(acc_muPlus))
  print("mu+ acceptance max:", np.max(acc_muPlus))
  print("mu+ acceptance mean:", np.mean(acc_muPlus))

  print("mu- acceptance min:", np.min(acc_muMinus))
  print("mu- acceptance max:", np.max(acc_muMinus))
  print("mu- acceptance mean:", np.mean(acc_muMinus))
  # Error mu+ and mu- 
  acc_err_muPlus, acc_err_muMinus = acc_error_3d(arrays_int_t, project="t") 
  # Make the plot 
  plot_acceptance_project_t(acc_muPlus, acc_muMinus, acc_err_muPlus, acc_err_muMinus)

# **********************************
def wrapper_project_phi():
  shape_int_phi = (n_t_bins, n_Q2_bins, n_nu_bins)
  arrays_int_phi = init_acceptance_arrays(shape_int_phi)
  fill_weights_3d(tree=tree_gen_muPlus, charge="muPlus", data_type="gen", arrays=arrays_int_phi, project="phi")
  fill_weights_3d(tree=tree_gen_muMinus, charge="muMinus", data_type="gen", arrays=arrays_int_phi, project="phi")
  fill_weights_3d(tree=tree_rec_muPlus, charge="muPlus", data_type="rec", arrays=arrays_int_phi, project="phi")
  fill_weights_3d(tree=tree_rec_muMinus, charge="muMinus", data_type="rec", arrays=arrays_int_phi, project="phi")
  # acceptance mu+ 
  acc_muPlus = np.zeros_like(arrays_int_phi["rec_plus"])
  mask_muPlus = arrays_int_phi["gen_plus"] != 0
  acc_muPlus[mask_muPlus] = ( arrays_int_phi["rec_plus"][mask_muPlus] / arrays_int_phi["gen_plus"][mask_muPlus]) 
  # acceptance mu- 
  acc_muMinus = np.zeros_like(arrays_int_phi["rec_minus"])
  mask_muMinus = arrays_int_phi["gen_minus"] != 0
  acc_muMinus[mask_muMinus] = ( arrays_int_phi["rec_minus"][mask_muMinus] / arrays_int_phi["gen_minus"][mask_muMinus]) 
  print("mu+ acceptance min:", np.min(acc_muPlus))
  print("mu+ acceptance max:", np.max(acc_muPlus))
  print("mu+ acceptance mean:", np.mean(acc_muPlus))
  print("mu- acceptance min:", np.min(acc_muMinus))
  print("mu- acceptance max:", np.max(acc_muMinus))
  print("mu- acceptance mean:", np.mean(acc_muMinus))
  # Error mu+ and mu- 
  acc_err_muPlus, acc_err_muMinus = acc_error_3d(arrays_int_phi, project="phi") 
  # Make the plot 
  plot_acceptance_project_phi(acc_muPlus, acc_muMinus, acc_err_muPlus, acc_err_muMinus)


# *******************************************************************
# *                       *** MAIN ***                              *
# *******************************************************************
# **********************************
# Main function - exclude or include functions here 
def main():
  # Make the default phase space plots with the loose binning 
  #default_plots()

  # 3D Acceptance - projected over |t| 
  #wrapper_project_t()
    
  # 3D Acceptance - projected over phi 
  # wrapper_project_phi()

  
  # Final 4D shape: [Q2][nu][t][phi]
  shape_4d = (
    len(Q2_edges_tight) - 1,
    len(nu_edges_tight) - 1,
    len(t_edges) - 1,
    len(phi_edges) - 1
  )
  arrays_4d = init_acceptance_arrays(shape_4d)
  
  # Fill the arrays 
  fill_weights_4d(tree_rec_muPlus, "muPlus", "rec", arrays_4d)
  fill_weights_4d(tree_gen_muPlus, "muPlus", "gen", arrays_4d)
  fill_weights_4d(tree_rec_muMinus, "muMinus", "rec", arrays_4d)
  fill_weights_4d(tree_gen_muMinus, "muMinus", "gen", arrays_4d)
  # Get the acceptance 
  # mu+ 
  acc_muPlus = np.zeros_like(arrays_4d["rec_plus"])
  mask_muPlus = arrays_4d["gen_plus"] != 0

  acc_muPlus[mask_muPlus] = (
    arrays_4d["rec_plus"][mask_muPlus] /
    arrays_4d["gen_plus"][mask_muPlus]
  )
  # mu-
  acc_muMinus = np.zeros_like(arrays_4d["rec_minus"])
  mask_muMinus = arrays_4d["gen_minus"] != 0

  acc_muMinus[mask_muMinus] = (
    arrays_4d["rec_minus"][mask_muMinus] /
    arrays_4d["gen_minus"][mask_muMinus]
  )
  print("mu+ acceptance min:", np.min(acc_muPlus))
  print("mu+ acceptance max:", np.max(acc_muPlus))
  print("mu+ acceptance mean:", np.mean(acc_muPlus))
  print("mu- acceptance min:", np.min(acc_muMinus))
  print("mu- acceptance max:", np.max(acc_muMinus))
  print("mu- acceptance mean:", np.mean(acc_muMinus))
  # Error mu+ and mu- 
  acc_err_muPlus, acc_err_muMinus = acc_error_4d(arrays_4d)

  with open("acc_G_P04.pkl", "wb") as f:
    pickle.dump(acc_muPlus, f)
    pickle.dump(acc_muMinus, f)
  print("Saved to acc_G_P04.pkl.")


  """
  A_global = (
    np.sum(arrays_4d["rec_plus"]) /
    np.sum(arrays_4d["gen_plus"])
  )
  print("GLOBAL acceptance:", A_global)

  # Plot by bin 
  plot_acceptance_by_tBin_4D(arrays_4d, acc_err_muPlus, acc_err_muMinus, 0)
  plot_acceptance_by_tBin_4D(arrays_4d, acc_err_muPlus, acc_err_muMinus, 1)
  plot_acceptance_by_tBin_4D(arrays_4d, acc_err_muPlus, acc_err_muMinus, 2)
  plot_acceptance_by_tBin_4D(arrays_4d, acc_err_muPlus, acc_err_muMinus, 3)
  """ 

if __name__ == "__main__":
  main()


