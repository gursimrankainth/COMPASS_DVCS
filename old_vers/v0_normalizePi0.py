import os
import glob
import ROOT
import numpy as np
from array import array
import datetime
from dvcs_constants_2016 import *

# **********************************
print("/!\ Make sure to use correct directories if calculating r_lepto")

# **********************************
# Real data 
real_dir = "/eos/user/g/gkainth/real/"
#real_dir = "/eos/user/g/gkainth/r_lepto/real/"
real_files = [os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root")]

# HEPGEN Pi0 MC Data (Reconstructed data)

hep_dir = "/eos/user/g/gkainth/HepgenPi0/"
#hep_dir = "/eos/user/g/gkainth/r_lepto/hepgen_pi0/"
hep_files = [os.path.join(hep_dir, "filtered_P04_muPlus.root"), os.path.join(hep_dir, "filtered_P04_muMinus.root"),
             os.path.join(hep_dir, "filtered_P05_muPlus.root"), os.path.join(hep_dir, "filtered_P05_muMinus.root"),
             os.path.join(hep_dir, "filtered_P06_muPlus.root"), os.path.join(hep_dir, "filtered_P06_muMinus.root"),
             os.path.join(hep_dir, "filtered_P07_muPlus.root"), os.path.join(hep_dir, "filtered_P07_muMinus.root"),
             os.path.join(hep_dir, "filtered_P08_muPlus.root"), os.path.join(hep_dir, "filtered_P08_muMinus.root"),
             os.path.join(hep_dir, "filtered_P09_muPlus.root"), os.path.join(hep_dir, "filtered_P09_muMinus.root")]

# LEPTO Visible Pi0 MC Data (Reconstructed data)
lepVis_dir = "/eos/user/g/gkainth/LeptoPi0/"
#lepVis_dir = "/eos/user/g/gkainth/r_lepto/lepto_pi0/"
lepVis_files = [os.path.join(lepVis_dir, "filtered_P04_muPlus.root"), os.path.join(lepVis_dir, "filtered_P04_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P05_muPlus.root"), os.path.join(lepVis_dir, "filtered_P05_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P06_muPlus.root"), os.path.join(lepVis_dir, "filtered_P06_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P07_muPlus.root"), os.path.join(lepVis_dir, "filtered_P07_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P08_muPlus.root"), os.path.join(lepVis_dir, "filtered_P08_muMinus.root"),
                os.path.join(lepVis_dir, "filtered_P09_muPlus.root"), os.path.join(lepVis_dir, "filtered_P09_muMinus.root")]

# LEPTO Exclusive Pi0 MC Data (Reconstructed data)
lepExcl_dir = "/eos/user/g/gkainth/old/LeptoExcl_new/"
#lepExcl_dir = "/eos/user/g/gkainth/r_lepto/lepExcl/"
lepExcl_files = [os.path.join(lepExcl_dir, "filtered_P04_muPlus.root"), os.path.join(lepExcl_dir, "filtered_P04_muMinus.root"),
                 os.path.join(lepExcl_dir, "filtered_P05_muPlus.root"), os.path.join(lepExcl_dir, "filtered_P05_muMinus.root"),
                 os.path.join(lepExcl_dir, "filtered_P06_muPlus.root"), os.path.join(lepExcl_dir, "filtered_P06_muMinus.root"),
                 os.path.join(lepExcl_dir, "filtered_P07_muPlus.root"), os.path.join(lepExcl_dir, "filtered_P07_muMinus.root"),
                 os.path.join(lepExcl_dir, "filtered_P08_muPlus.root"), os.path.join(lepExcl_dir, "filtered_P08_muMinus.root"),
                 os.path.join(lepExcl_dir, "filtered_P09_muPlus.root"), os.path.join(lepExcl_dir, "filtered_P09_muMinus.root")]


# **********************************
# ***     Exclusivity Plots      ***
# **********************************
# /!\ PLOTS ARE MADE WITH ALL DATA, NO GENERAL NORMALIZATION FACTORS APPLIED YET
# Real data exclusivity variable plots  
def real_excl_plots(nBins=10):
  # Create histograms
  h_real_pt = ROOT.TH1D("h_real_pt", "Real Data; #Delta p_{T} [GeV/c]; Counts", nBins, -0.3, 0.3)
  h_real_phi = ROOT.TH1D("h_real_phi", "Real Data; #Delta #phi [rad]; Counts", nBins, -0.4, 0.4)
  h_real_Z = ROOT.TH1D("h_real_Z", "Real Data; #Delta Z [cm]; Counts", nBins, -15, 15)
  h_real_mm = ROOT.TH1D("h_real_mm", "Real Data; #Delta M_{miss}^{2} [GeV^{2}/c^{4}]; Counts", nBins, -0.3, 0.3)
  h_inv_mass = ROOT.TH1D("h_inv_mass", "Real Data; M_{#gamma#gamma}; Counts", nBins, 0, 0.3)

  # Create TChain and add files
  chain = ROOT.TChain("USR970_pi0")
  for f in real_files:
    chain.Add(f)

  # Loop over events
  for event in chain:
    M_p = 0.938272  # proton mass
    targ_TL = ROOT.TLorentzVector(0, 0, 0, M_p)

    # Reconstruct pi0 4-vector
    pi0_TL = event.gamma_TL + event.gammaLow_TL
    h_inv_mass.Fill(pi0_TL.M())

    # Recalculate proton spec vector for pi0 topology
    p_spec_TL = targ_TL + event.inMu_TL - event.outMu_TL - pi0_TL

    # Δpₜ 
    delta_pt = event.p_camera_TL.Pt() - p_spec_TL.Pt()
    h_real_pt.Fill(delta_pt)

    # Δφ: difference in azimuthal angle of proton (reconstructed vs spec)
    delta_phi = event.p_camera_TL.Phi() - p_spec_TL.Phi()
    # Wrap into [-pi, pi]
    while delta_phi > ROOT.TMath.Pi():
      delta_phi -= 2*ROOT.TMath.Pi()
    while delta_phi < -ROOT.TMath.Pi():
      delta_phi += 2*ROOT.TMath.Pi()
    h_real_phi.Fill(delta_phi)

    # ΔZ: difference in z position in ring A (measured vs. interpolated)
    h_real_Z.Fill(event.delta_Z)

    # Missing mass squared: (target + inMu - outMu - pi0)^2
    delta_M2x = (event.p_camera_TL - p_spec_TL) * (event.p_camera_TL - p_spec_TL);
    h_real_mm.Fill(delta_M2x)

  # Create canvas and divide into 2x2 pads
  c = ROOT.TCanvas("c_real_excl", "Exclusivity variables", 1200, 1000)
  c.Divide(2, 2)

  # Draw histograms in subpads
  c.cd(1)
  h_real_pt.Draw()
  c.cd(2)
  h_real_phi.Draw()
  c.cd(3)
  h_real_Z.Draw()
  c.cd(4)
  h_real_mm.Draw()

  # Update and save
  c.Update()
  c.SaveAs("real_visPi0_exclusivity.png")

  c1 = ROOT.TCanvas("c_inv_mass", "Invairant Mass", 800, 600)
  h_inv_mass.Draw()
  c1.Update()
  c1.SaveAs("real_inv_mass.png")

  return (h_real_pt, h_real_phi, h_real_Z, h_real_mm)

# **********************************
# Exclusive pi0 (hepgen) data exclusivity variable plots 
def hep_excl_plots(nBins=10):
  # Create histograms
  h_hep_pt = ROOT.TH1D("h_hep_pt", "Exclusive #pi^{0} Data; #Delta p_{T} [GeV/c]; Counts", nBins, -0.3, 0.3)
  h_hep_phi = ROOT.TH1D("h_hep_phi", "Exclusive #pi^{0} Data; #Delta #phi [rad]; Counts", nBins, -0.4, 0.4)
  h_hep_Z = ROOT.TH1D("h_hep_Z", "Exclusive #pi^{0} Data; #Delta Z [cm]; Counts", nBins, -15, 15)
  h_hep_mm = ROOT.TH1D("h_hep_mm", "Exclusive #pi^{0} Data; #Delta M_{miss}^{2} [GeV^{2}/c^{4}]; Counts", nBins, -0.3, 0.3)

  # Create TChain and add files
  chain = ROOT.TChain("USR970_pi0")
  for f in hep_files:
    chain.Add(f)

  # Loop over events
  for event in chain:
    M_p = 0.938272  # proton mass
    targ_TL = ROOT.TLorentzVector(0, 0, 0, M_p)
    weight = event.weight_all

    # Reconstruct pi0 4-vector
    pi0_TL = event.gamma_TL + event.gammaLow_TL

    # Recalculate proton spec vector for pi0 topology
    p_spec_TL = targ_TL + event.inMu_TL - event.outMu_TL - pi0_TL

    # Δpₜ 
    delta_pt = event.p_camera_TL.Pt() - p_spec_TL.Pt()
    h_hep_pt.Fill(delta_pt, weight)

    # Δφ: difference in azimuthal angle of proton (reconstructed vs spec)
    delta_phi = event.p_camera_TL.Phi() - p_spec_TL.Phi()
    # Wrap into [-pi, pi]
    while delta_phi > ROOT.TMath.Pi():
      delta_phi -= 2*ROOT.TMath.Pi()
    while delta_phi < -ROOT.TMath.Pi():
      delta_phi += 2*ROOT.TMath.Pi()
    h_hep_phi.Fill(delta_phi, weight)

    # ΔZ: difference in z position in ring A (measured vs. interpolated)
    h_hep_Z.Fill(event.delta_Z, weight)

    # Missing mass squared: (target + inMu - outMu - pi0)^2
    delta_M2x = (event.p_camera_TL - p_spec_TL) * (event.p_camera_TL - p_spec_TL);
    h_hep_mm.Fill(delta_M2x, weight)

  # Create canvas and divide into 2x2 pads
  c = ROOT.TCanvas("c_hep_excl", "Exclusivity variables", 1200, 1000)
  c.Divide(2, 2)

  # Draw histograms in subpads
  c.cd(1)
  h_hep_pt.Draw("hist")
  c.cd(2)
  h_hep_phi.Draw("hist")
  c.cd(3)
  h_hep_Z.Draw("hist")
  c.cd(4)
  h_hep_mm.Draw("hist")

  # Update and save
  c.Update()
  c.SaveAs("hep_visPi0_exclusivity.png")

  return (h_hep_pt, h_hep_phi, h_hep_Z, h_hep_mm)

# **********************************
# Inclusive pi0 (lepto) data exclusivity variable plots 
def lep_excl_plots(nBins=10):
  # Create histograms
  h_lep_pt = ROOT.TH1D("h_lep_pt", "Inclusive #pi^{0} Data; #Delta p_{T} [GeV/c]; Counts", nBins, -0.3, 0.3)
  h_lep_phi = ROOT.TH1D("h_lep_phi", "Inclusive #pi^{0} Data; #Delta #phi [rad]; Counts", nBins, -0.4, 0.4)
  h_lep_Z = ROOT.TH1D("h_lep_Z", "Inclusive #pi^{0} Data; #Delta Z [cm]; Counts", nBins, -15, 15)
  h_lep_mm = ROOT.TH1D("h_lep_mm", "Inclusive #pi^{0} Data; #Delta M_{miss}^{2} [GeV^{2}/c^{4}]; Counts", nBins, -0.3, 0.3)

  # Create TChain and add files
  chain = ROOT.TChain("USR970_pi0")
  for f in lepExcl_files:
    chain.Add(f)

  # Loop over events
  for event in chain:
    # Δpₜ 
    h_lep_pt.Fill(event.delta_pt)

    # Δφ: difference in azimuthal angle of proton (reconstructed vs spec)
    h_lep_phi.Fill(event.delta_phi)

    # ΔZ: difference in z position in ring A (measured vs. interpolated)
    h_lep_Z.Fill(event.delta_Z)

    # Missing mass squared: (target + inMu - outMu - pi0)^2
    h_lep_mm.Fill(event.M2x)

  # Create canvas and divide into 2x2 pads
  c = ROOT.TCanvas("c_lep_excl", "Exclusivity variables", 1200, 1000)
  c.Divide(2, 2)

  # Draw histograms in subpads
  c.cd(1)
  h_lep_pt.Draw()
  c.cd(2)
  h_lep_phi.Draw()
  c.cd(3)
  h_lep_Z.Draw()
  c.cd(4)
  h_lep_mm.Draw()

  # Update and save
  c.Update()
  c.SaveAs("lep_exclPi0_exclusivity.png")

  return (h_lep_pt, h_lep_phi, h_lep_Z, h_lep_mm)


# **********************************
# ***     c_HEP and c_LEP        ***
# **********************************
# Get the general scale factors (one per beam charge for each mc sample)
def gen_scaleFacs(periods, debug_mode=False):
  hep_muPlus_lst = [] 
  hep_muMinus_lst = [] 
  lep_muPlus_lst = []
  lep_muMinus_lst = [] 

  def count_events(period, file_list, debug_mode=debug_mode, use_weight=False, use_filename_charge=False):
    chain = ROOT.TChain("USR970_pi0")

    count_muPlus = 0.0
    count_muMinus = 0.0

    for f in file_list:
      if f"_{period}" in f:
        if debug_mode: 
          print(f)
        chain.Add(f)

        # Determine charge from filename if requested
        file_charge = None
        if use_filename_charge:
          if "muPlus" in f:
            file_charge = +1
          elif "muMinus" in f:
            file_charge = -1

        for event in chain:
          w = event.weight_all if use_weight else 1.0

          if use_filename_charge:
            if file_charge == +1:
              count_muPlus += w
            elif file_charge == -1:
              count_muMinus += w
          else:
            if event.Q_beam == 1:
              count_muPlus += w
            elif event.Q_beam == -1:
              count_muMinus += w

        chain.Reset()  # important: avoid re-looping previous files

    return count_muPlus, count_muMinus

  for period in periods: 
    if debug_mode:
      print(period)
      
    real_muPlus, real_muMinus = count_events(period, real_files, use_weight=False)
    hep_muPlus, hep_muMinus = count_events(period, hep_files, use_weight=True, use_filename_charge=True)
    lep_muPlus, lep_muMinus = count_events(period, lepVis_files, use_weight=False, use_filename_charge=True)

    # Normalization constants
    c_hep_muPlus  = real_muPlus / hep_muPlus  if hep_muPlus  != 0 else 0.0
    hep_muPlus_lst.append(c_hep_muPlus)
    c_hep_muMinus = real_muMinus / hep_muMinus if hep_muMinus != 0 else 0.0
    hep_muMinus_lst.append(c_hep_muMinus)
    c_lep_muPlus  = real_muPlus / lep_muPlus  if lep_muPlus  != 0 else 0.0
    lep_muPlus_lst.append(c_lep_muPlus)
    c_lep_muMinus = real_muMinus / lep_muMinus if lep_muMinus != 0 else 0.0
    lep_muMinus_lst.append(c_lep_muMinus)

    if debug_mode:
      print("HEPGEN:", c_hep_muPlus, c_hep_muMinus)
      print("LEPTO :", c_lep_muPlus, c_lep_muMinus)

  return hep_muPlus_lst, lep_muPlus_lst, hep_muMinus_lst, lep_muMinus_lst

# **********************************
#  Write calculated values to Python file for storage
def write_constants_pi0(const_hep_mup, const_lep_mup, const_hep_mum, const_lep_mum, r_lepto, year=2016):
  constants_file = f"/afs/cern.ch/user/g/gkainth/dvcs_constants_{year}.py" 
  if not os.path.exists(constants_file):
    print(f"No constants file found for {year}: {constants_file}")
    return

  source_file = os.path.basename(__file__)
  today = datetime.date.today().isoformat()  # "YYYY-MM-DD"

  with open(constants_file, "r") as f:
    lines = f.readlines()
  
  # Safety checks
  assert "CBH_MUPLUS" in lines[8]
  assert "CBH_MUMINUS" in lines[9]

  lines[11] = f'SOURCE_FILE = "{source_file} [{today}]"\n'
  lines[12] = f"CPI0_HEP_MUPLUS = {const_hep_mup}\n"
  lines[13] = f"CPI0_HEP_MUMINUS = {const_hep_mum}\n"
  lines[14] = f"CPI0_LEP_MUPLUS = {const_lep_mup}\n"
  lines[15] = f"CPI0_LEP_MUMINUS = {const_lep_mum}\n"
  lines[16] = f"R_LEPTO = {r_lepto}\n"

  with open(constants_file, "w") as f:
    f.writelines(lines)

  print("Values written to file:", constants_file)


# **********************************
# ***         r_LEPTO            ***
# **********************************
def get_r_lepto(real_hists, hep_hists, lep_hists, excl_var="delta_pt"):
  if excl_var == "delta_pt":
    h_real = real_hists[0]
    h_hep = hep_hists[0]
    h_lep = lep_hists[0]
  elif excl_var == "delta_phi":
    h_real = real_hists[1]
    h_hep = hep_hists[1]
    h_lep = lep_hists[1]
  elif excl_var == "delta_Z":
    h_real = real_hists[2]
    h_hep = hep_hists[2]
    h_lep = lep_hists[2]
  elif excl_var == "delta_M2x":
    h_real = real_hists[3]
    h_hep = hep_hists[3]
    h_lep = lep_hists[3]
  
  h_real.SetDirectory(0)
  h_hep.SetDirectory(0)
  h_lep.SetDirectory(0)

  h_real.Sumw2()
  h_hep.Sumw2()
  h_lep.Sumw2()

  # Prepare MC array
  mc_array = ROOT.TObjArray(2)
  mc_array.Add(h_hep)  # 0 → HEPGEN
  mc_array.Add(h_lep)  # 1 → LEPTO

  fitter = ROOT.TFractionFitter(h_real, mc_array)
  fitter.Constrain(1, 0.0, 1.0)  # LEPTO fraction [0,1]

  status = fitter.Fit()
  if status != 0:
      print(f"WARNING: TFractionFitter did not converge, status = {status}")

  # Get fitted fractions
  R_lepto   = array('d', [0.0])
  err_lepto = array('d', [0.0])
  R_hep     = array('d', [0.0])
  err_hep   = array('d', [0.0])

  fitter.GetResult(1, R_lepto, err_lepto)
  fitter.GetResult(0, R_hep, err_hep)

  fit_hist = ROOT.TH1D(fitter.GetPlot())
  fit_hist.SetDirectory(0)

  return R_lepto[0], err_lepto[0], R_hep[0], err_hep[0], status, fit_hist

# **********************************
# Plot with fit results (mc only)
def compare_mc_only(fit_hist, excl_var="delta_pt", r_lepto=0.4, nBins=16):

  # --- Axis setup ---
  if excl_var == "delta_pt":
    xmin, xmax = -0.3, 0.3
    xTitle = "#Delta p_{T} [GeV/c]"
  elif excl_var == "delta_phi":
    xmin, xmax = -0.4, 0.4
    xTitle = "#Delta #phi [rad]"
  elif excl_var == "delta_Z":
    xmin, xmax = -15, 15
    xTitle = "#Delta Z_{A} [cm]"
  elif excl_var == "delta_M2x":
    xmin, xmax = -0.3, 0.3
    xTitle = "#Delta M_{miss}^{2} [GeV^{2}/c^{4}]"

  # --- Histograms ---
  h_hep  = ROOT.TH1D("h_hep", f"Exclusive #pi0 MC; {xTitle}; Counts", nBins, xmin, xmax)
  h_lep  = ROOT.TH1D("h_lep", f"Inclusive #pi0 MC; {xTitle}; Counts", nBins, xmin, xmax)

  h_hep.Sumw2()
  h_lep.Sumw2()

  M_p = 0.938272
  targ_TL = ROOT.TLorentzVector(0, 0, 0, M_p)

  # --- Fill HEPGEN ---
  chain_hep = ROOT.TChain("USR970_pi0")
  for f in hep_files:
    chain_hep.Add(f)

  for event in chain_hep:
    weight = event.weight_all
    pi0_TL = event.gamma_TL + event.gammaLow_TL
    p_spec_TL = targ_TL + event.inMu_TL - event.outMu_TL - pi0_TL

    if excl_var == "delta_pt":
      val = event.p_camera_TL.Pt() - p_spec_TL.Pt()
    elif excl_var == "delta_phi":
      val = event.p_camera_TL.Phi() - p_spec_TL.Phi()
      while val > ROOT.TMath.Pi(): val -= 2*ROOT.TMath.Pi()
      while val < -ROOT.TMath.Pi(): val += 2*ROOT.TMath.Pi()
    elif excl_var == "delta_Z":
      val = event.delta_Z
    elif excl_var == "delta_M2x":
      val = (event.p_camera_TL - p_spec_TL) * (event.p_camera_TL - p_spec_TL)

    h_hep.Fill(val, weight)

  # --- Fill LEPTO ---
  chain_lep = ROOT.TChain("USR970_pi0")
  for f in lepExcl_files:
    chain_lep.Add(f)

  for event in chain_lep:
    weight = getattr(event, "weight_all", 1.0)
    if excl_var == "delta_pt":
      val = event.delta_pt
    elif excl_var == "delta_phi":
      val = event.delta_phi
    elif excl_var == "delta_Z":
      val = event.delta_Z
    elif excl_var == "delta_M2x":
      val = event.M2x
    h_lep.Fill(val, weight)

  # --- Normalize shapes to unity ---
  h_hep_scaled = h_hep.Clone("h_hep_scaled")
  h_lep_scaled = h_lep.Clone("h_lep_scaled")
  if h_hep_scaled.Integral() > 0: h_hep_scaled.Scale(1.0 / h_hep_scaled.Integral())
  if h_lep_scaled.Integral() > 0: h_lep_scaled.Scale(1.0 / h_lep_scaled.Integral())

  # --- Apply fractions once and build total mixture ---
  h_hep_scaled.Scale(1 - r_lepto)
  h_lep_scaled.Scale(r_lepto)

  h_total = h_hep_scaled.Clone("h_total")
  h_total.Add(h_lep_scaled)

  # --- Scale total mixture to match fit_hist integral ---
  if fit_hist.Integral() > 0 and h_total.Integral() > 0:
    scale_factor = fit_hist.Integral() / h_total.Integral()
    #h_total.Scale(scale_factor)
    h_hep_scaled.Scale(scale_factor)
    h_lep_scaled.Scale(scale_factor)

  # --- Create fit points with Poisson errors ---
  graph_fit = ROOT.TGraphErrors()
  for i in range(1, fit_hist.GetNbinsX() + 1):
    x = fit_hist.GetBinCenter(i)
    y = fit_hist.GetBinContent(i)
    ex = 0.0
    ey = np.sqrt(y) if y > 0 else 0.0
    graph_fit.SetPoint(i-1, x, y)
    graph_fit.SetPointError(i-1, ex, ey)

  # --- Axis scaling ---
  max_val = max(h_total.GetMaximum(), fit_hist.GetMaximum())
  h_total.SetMaximum(1.2 * max_val)
  h_total.SetMinimum(0)

  # --- Canvas and drawing ---
  c = ROOT.TCanvas("c_mc_only", "MC comparison", 800, 600)

  h_total.SetLineColor(ROOT.kBlack)
  h_total.SetLineWidth(2)
  h_total.Draw("hist")

  h_hep_scaled.SetLineColor(ROOT.kBlue)
  h_hep_scaled.SetLineWidth(2)
  h_hep_scaled.SetFillColorAlpha(ROOT.kBlue, 0.15)
  h_hep_scaled.Draw("hist same")

  h_lep_scaled.SetLineColor(ROOT.kRed)
  h_lep_scaled.SetLineWidth(2)
  h_lep_scaled.SetFillColorAlpha(ROOT.kRed, 0.15)
  h_lep_scaled.Draw("hist same")

  graph_fit.SetMarkerStyle(20)
  graph_fit.SetMarkerSize(0.9)
  graph_fit.SetMarkerColor(ROOT.kViolet-5)
  graph_fit.SetLineColor(ROOT.kViolet-5)
  graph_fit.SetLineWidth(1)
  graph_fit.Draw("P same")

  # --- Legend ---
  leg = ROOT.TLegend(0.15, 0.72, 0.40, 0.85)
  leg.AddEntry(h_hep_scaled, "HEPGEN MC", "f")
  leg.AddEntry(h_lep_scaled, "LEPTO MC", "f")
  leg.AddEntry(h_total, "MC mixture (scaled to fit)", "l")
  leg.AddEntry(graph_fit, "Fit", "p")
  leg.Draw()

  c.Update()
  c.SaveAs(f"mc_only_{excl_var}.png")

  return h_hep_scaled, h_lep_scaled, h_total



"""real_hists = real_excl_plots(nBins=16)
hep_hists = hep_excl_plots(nBins=16)
lep_hists = lep_excl_plots(nBins=16)

r_lepto, err_lepto, r_hep, err_hep, status, fit_hist = get_r_lepto(real_hists, hep_hists, lep_hists, excl_var="delta_M2x")
print("delta_M2x:", r_lepto, err_lepto, r_hep, err_hep)
compare_mc_only(fit_hist, r_lepto=r_lepto, excl_var="delta_M2x") 

 r_lepto, err_lepto, r_hep, err_hep, status, fit_hist = get_r_lepto(real_hists, hep_hists, lep_hists, excl_var="delta_phi")
print("delta_phi:", r_lepto, err_lepto, r_hep, err_hep)
compare_mc_only(fit_hist, r_lepto=r_lepto, excl_var="delta_phi")

r_lepto, err_lepto, r_hep, err_hep, status, fit_hist = get_r_lepto(real_hists, hep_hists, lep_hists, excl_var="delta_pt")
print("delta_pt:", r_lepto, err_lepto, r_hep, err_hep)
compare_mc_only(fit_hist, r_lepto=r_lepto, excl_var="delta_pt")

r_lepto, err_lepto, r_hep, err_hep, status, fit_hist = get_r_lepto(real_hists, hep_hists, lep_hists, excl_var="delta_Z")
print("delta_Z:", r_lepto, err_lepto, r_hep, err_hep)
compare_mc_only(fit_hist, r_lepto=r_lepto, excl_var="delta_Z")"""


hep_muPlus_lst, lep_muPlus_lst, hep_muMinus_lst, lep_muMinus_lst = gen_scaleFacs(PERIODS, debug_mode=False)
r_lepto = 0.4 
write_constants_pi0(const_hep_mup=hep_muPlus_lst, const_hep_mum=hep_muMinus_lst,
                    const_lep_mup=lep_muPlus_lst, const_lep_mum=lep_muMinus_lst,
                    r_lepto=r_lepto)