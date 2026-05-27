import os
import glob
import ROOT
import numpy as np
from array import array
import datetime

# **********************************
# HEPGEN BH MC Data (generated data)
gen_dir = "/Users/gursimran/cern/2016_data/BH/"
mcGen_files = [[os.path.join(gen_dir, "gen_P04_muPlus.root"), os.path.join(gen_dir, "gen_P04_muMinus.root")],
               [os.path.join(gen_dir, "gen_P05_muPlus.root"), os.path.join(gen_dir, "gen_P05_muMinus.root")],
               [os.path.join(gen_dir, "gen_P06_muPlus.root"), os.path.join(gen_dir, "gen_P06_muMinus.root")],
               [os.path.join(gen_dir, "gen_P07_muPlus.root"), os.path.join(gen_dir, "gen_P07_muMinus.root")],
               [os.path.join(gen_dir, "gen_P08_muPlus.root"), os.path.join(gen_dir, "gen_P08_muMinus.root")],
               [os.path.join(gen_dir, "gen_P09_muPlus.root"), os.path.join(gen_dir, "gen_P09_muMinus.root")]]

# HEPGEN BH MC Data (reconstructed data)
rec_dir = "/Users/gursimran/cern/2016_data/BH/"
mcRec_files = [[os.path.join(rec_dir, "filtered_P04_muPlus.root"), os.path.join(rec_dir, "filtered_P04_muMinus.root")],
               [os.path.join(rec_dir, "filtered_P05_muPlus.root"), os.path.join(rec_dir, "filtered_P05_muMinus.root")],
               [os.path.join(rec_dir, "filtered_P06_muPlus.root"), os.path.join(rec_dir, "filtered_P06_muMinus.root")],
               [os.path.join(rec_dir, "filtered_P07_muPlus.root"), os.path.join(rec_dir, "filtered_P07_muMinus.root")],
               [os.path.join(rec_dir, "filtered_P08_muPlus.root"), os.path.join(rec_dir, "filtered_P08_muMinus.root")],
               [os.path.join(rec_dir, "filtered_P09_muPlus.root"), os.path.join(rec_dir, "filtered_P09_muMinus.root")]]

# Real data 
real_dir = "/Users/gursimran/cern/2016_data/real/"
real_files = [[os.path.join(real_dir, "filtered_P04.root")],
              [os.path.join(real_dir, "filtered_P05.root")],
              [os.path.join(real_dir, "filtered_P06.root")],
              [os.path.join(real_dir, "filtered_P07.root")],
              [os.path.join(real_dir, "filtered_P08.root")],
              [os.path.join(real_dir, "filtered_P09.root")]]

# Flux files 
flux_dir = "/Users/gursimran/cern/configs/flux_Johannes/2016/flux_files"
#flux_dir = "/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files_slot7.1"


# **********************************
# Integrated beam luminosity for each period (Always use RT corrected values)
period_run_ranges = {
  "P04": (272340, 273116),
  "P05": (273125, 273605),
  "P06": (273666, 274200),
  "P07": (274508, 274901),
  "P08": (274946, 275393),
  "P09": (275478, 275908),
}

def load_bad_spills(period, bad_spill_dir="/Users/gursimran/cern/configs/bad_spill/"):
  # Find matching bad spill file for the given period
  pattern = os.path.join(bad_spill_dir, f"P0{period}*bad_spill.lst")
  print(pattern)
  matching_files = glob.glob(pattern)
  if len(matching_files) == 0:
    print("Bad spill file not found")
    return set()

  bad_spill_file = matching_files[0]
  #print(f"Using bad spill file: {bad_spill_file} ... ")

  bad_spills = set()

  # Read bad spill file
  with open(bad_spill_file, "r") as f:
    for line in f:
      parts = line.strip().split()
      if len(parts) < 9:
        continue
      run, spill = int(parts[0]), int(parts[1])
      flux, LAST, LT, MT, OT, RICH, ECAL, empty = map(int, parts[2:])
      if flux == 1 or LT == 1 or MT == 1 or OT == 1 or ECAL == 1 or empty == 1:
        bad_spills.add((run, spill))

  #print(f"Loaded {len(bad_spills)} bad spills ...")
  return bad_spills

def get_flux(flux_dir=flux_dir, period=9):
  period_name = "P0" + str(period)
  run_min, run_max = period_run_ranges[period_name]

  # Load bad spills
  bad_spills = load_bad_spills(period)
  #print(f"Loaded {len(bad_spills)} bad spills ...")

  flux_muPlus = 0.0
  flux_muMinus = 0.0
  spill_count = 0

  # Loop over flux files
  for fname in sorted(os.listdir(flux_dir)):
    if not fname.startswith("flux_run") or not fname.endswith(".txt"):
      continue

    run_number = int(fname.replace("flux_run", "").replace(".txt", ""))
    if not (run_min <= run_number <= run_max):
      continue

    # Read flux file
    with open(os.path.join(flux_dir, fname), 'r') as f:
      for line in f.readlines()[2:]:
        entries = line.split()
        if len(entries) < 11:
          continue

        spill = int(entries[0])
        if (run_number, spill) in bad_spills:
          continue

        charge = int(entries[10])
        if charge == 0:
          continue

        spill_flux = float(entries[5])  # RT_flux(VDTcorr)
        spill_count += 1

        if charge == 1:
          flux_muPlus += spill_flux
        elif charge == -1:
          flux_muMinus += spill_flux

  return flux_muPlus, flux_muMinus, spill_count

def get_luminosity(intFlux_muPlus, intFlux_muMinus, period=9): 
  M_p = 1.00727646627       # g/mol molar proton mass
  rho_H2 = 0.070146         # g/cm3 proton density in liq. hydrogen
  Na = 6.02214076e23        # 1/mol Avogadro constant
  l = 240.0                 # cm target length

  idx = period - 4
  flux_p = float(intFlux_muPlus[idx])
  flux_m = float(intFlux_muMinus[idx])

  lumin_p = (rho_H2 * Na * l * flux_p) / M_p
  lumin_m = (rho_H2 * Na * l * flux_m) / M_p
  return lumin_p, lumin_m


# **********************************
# Integrated MC beam luminosity for each period 
def get_MC_luminosity(mc_files, period=9, charge="muMinus"):
  xSection_DVCS = 9.93938024e-34
  idx = period - 4

  file_name = mc_files[idx][0] if charge=="muPlus" else mc_files[idx][1]
  f = ROOT.TFile.Open(file_name)
  if not f or f.IsZombie():
    print(f"Could not open {file_name}")
    return

  tree = f.Get("USR970_gen")
  if not tree:
    print(f"Tree USR970_gen not found in {file_name}")
    return

  sum_weight = 0.0
  for i in range(tree.GetEntries()):
    tree.GetEntry(i)

    weight_DVCS = float(tree.weight_DVCS)
    sum_weight += weight_DVCS

  lumin = sum_weight / xSection_DVCS
  return lumin


# **********************************
# Plot the normalized MC BH data with the real data 
def plot_variable_allPeriods(
  const_muPlus,
  const_muMinus,
  var="nu_fit",
  x_min=0.0,
  x_max=150.0,
  n_bins=30,
  logy=False,
  logx=False,
  region="full",  # full, reference, interference, extraction
  scale_factor=None # cosmetic scale (applies to both)
  ):

  # *****************
  # Region cuts (based ONLY on nu)
  def in_region(tree):
    nu = tree.nu_fit

    if region == "full":
      return (10 < nu < 144)

    if region == "reference":       # 80–144
      return (80 < nu < 144)

    if region == "interference":    # 32–80
      return (32 < nu < 80)

    if region == "extraction":      # 10–32
      return (10 < nu < 32)

    # fallback (shouldn't happen)
    return True

  # *****************
  # Combined histograms
  hist_mc_all   = ROOT.TH1F(f"hist_mc_{var}",   f"MC Combined {var}",   n_bins, x_min, x_max)
  hist_data_all = ROOT.TH1F(f"hist_data_{var}", f"Data Combined {var}", n_bins, x_min, x_max)

  # *****************
  # Loop over all periods
  for period_idx in range(len(real_files)):

    norm_plus  = const_muPlus[period_idx]
    norm_minus = const_muMinus[period_idx]

    # ********** MC μ+ **********
    f_mc_plus = ROOT.TFile.Open(mcRec_files[period_idx][0])
    tree = f_mc_plus.Get("USR970_filtered")

    for i in range(tree.GetEntries()):
      tree.GetEntry(i)
      if in_region(tree):
        #weight_BH = getattr(tree, "weight_BH")
        weight_BH = getattr(tree, "weight_PAMBH")
        value = getattr(tree, var)
        if var == "t_fit":
          value = abs(value)
        hist_mc_all.Fill(value, norm_plus * weight_BH)
    f_mc_plus.Close()

    # ********** MC μ− **********
    f_mc_minus = ROOT.TFile.Open(mcRec_files[period_idx][1])
    tree = f_mc_minus.Get("USR970_filtered")

    for i in range(tree.GetEntries()):
      tree.GetEntry(i)
      if in_region(tree):
        #weight_BH = getattr(tree, "weight_BH")
        weight_BH = getattr(tree, "weight_PAMBH")
        value = getattr(tree, var)
        if var == "t_fit":
          value = abs(value)
        hist_mc_all.Fill(value, norm_minus * weight_BH)
    f_mc_minus.Close()

    # ********** Real data **********
    f_data = ROOT.TFile.Open(real_files[period_idx][0])
    tree = f_data.Get("USR970_filtered")

    for i in range(tree.GetEntries()):
      tree.GetEntry(i)
      if in_region(tree):
        value = getattr(tree, var)
        if var == "t_fit":
          value = abs(value)
        hist_data_all.Fill(value, 1.0)
    f_data.Close()

  # *****************
  # Style
  hist_data_all.SetLineColor(ROOT.kAzure+1)
  hist_data_all.SetMarkerColor(ROOT.kAzure+1)
  hist_data_all.SetMarkerStyle(20)
  hist_data_all.SetMarkerSize(0.7)

  hist_mc_all.SetLineColor(ROOT.kBlack)
  hist_mc_all.SetLineWidth(2)

  # *****************
  # Canvas
  c = ROOT.TCanvas(f"c_{var}_{region}", f"{var} ({region})", 800, 600)
  if logy:
    c.SetLogy()
  if logx:
    c.SetLogx()

  hist_mc_all.GetXaxis().SetTitle(var)

  # Build Y-axis label based on cosmetic scale
  y_label = "Events"
  if isinstance(scale_factor, (float, int)) and scale_factor != 1.0:
    if var == "nu_fit":
      y_label = f"Events / {scale_factor:.2f} GeV"
    elif var == "Q2_fit":
      y_label = f"Events / {scale_factor:.2f} (GeV/c)^2"
    elif var == "t_fit":
      y_label = f"Events / {scale_factor:.2f} GeV^2"
    else: 
      y_label = f"Events"

  hist_mc_all.GetYaxis().SetTitle(y_label)
  hist_data_all.GetYaxis().SetTitle(y_label)

  # Cosmetic scale factor applies to BOTH histograms
  if isinstance(scale_factor, (float, int)):
    print(f"Applying cosmetic scale factor: {scale_factor}")
    hist_mc_all.Scale(scale_factor)
    hist_data_all.Scale(scale_factor)
  
  # Scale for the y_axis
  if var == "Q2_fit":
    hist_mc_all.GetYaxis().SetRangeUser(1, hist_mc_all.GetMaximum()*2)
  else:
    hist_mc_all.GetYaxis().SetRangeUser(1, hist_mc_all.GetMaximum()*1.2)

  hist_mc_all.Draw("HIST")
  hist_data_all.Draw("PE SAME")

  # *****************
  # Legend (upper-left, semi-transparent)
  legend = ROOT.TLegend(0.15, 0.75, 0.35, 0.88)
  legend.SetBorderSize(0)
  legend.SetFillColorAlpha(0, 0.3)  # semi-transparent box
  legend.SetTextSize(0.035)

  legend.AddEntry(hist_data_all, "Data (all periods)", "pe")
  legend.AddEntry(hist_mc_all,   "MC (all periods)",   "l")

  #legend.Draw()

  # Save figure
  output_name = f"{var}_{region}_allPeriods.png"
  c.SaveAs(output_name)
  c.Close()

  print(f"Created: {output_name}")


# **********************************
# Write calculated values to Python file for storage
def write_constants_BH(intLum_muPlus, intLum_muMinus, const_muPlus, const_muMinus, year=2016):
  constants_file = f"/afs/cern.ch/user/g/gkainth/dvcs_constants_{year}.py" 
  if not os.path.exists(constants_file):
    print(f"No constants file found for {year}: {constants_file}")
    return

  source_file = os.path.basename(__file__)
  today = datetime.date.today().isoformat()  # "YYYY-MM-DD"

  with open(constants_file, "r") as f:
    lines = f.readlines()
  
  # Safety checks
  assert "SOURCE_FILE" in lines[5]
  assert "LUMINOSITY_MUPLUS" in lines[6]
  assert "LUMINOSITY_MUMINUS" in lines[7]

  lines[5] = f'SOURCE_FILE = "{source_file} [{today}]"\n'
  lines[6] = f"LUMINOSITY_MUPLUS = {intLum_muPlus}\n"
  lines[7] = f"LUMINOSITY_MUMINUS = {intLum_muMinus}\n"
  lines[8] = f"CBH_MUPLUS = {const_muPlus}\n"
  lines[9] = f"CBH_MUMINUS = {const_muMinus}\n"

  with open(constants_file, "w") as f:
    f.writelines(lines)

  print("Values written to file:", constants_file)


# **********************************
# Main function - exclude or include functions here 
def main(): 
  # *** Real data ***
  # flux 
  intFlux_muPlus = []
  intFlux_muMinus = []
  for i in range(4,10):
    flux_p, flux_m, n_spills = get_flux(period=i)
    intFlux_muPlus.append(flux_p)
    intFlux_muMinus.append(flux_m)

  # luminosity 
  intLum_muPlus = []
  intLum_muMinus = [] 
  for i in range(4,10):
    lumin_p, lumin_m = get_luminosity(intFlux_muPlus, intFlux_muMinus, period=i)
    print("Period:", i, ", Lumin_p:", lumin_p)
    print("Period:", i, ", Lumin_m:", lumin_m)
    intLum_muPlus.append(lumin_p)
    intLum_muMinus.append(lumin_m)

  # *** MC data ***
  # luminosity (using generated data) 
  intLumMC_muPlus = [] 
  intLumMC_muMinus = [] 
  for i in range(4,10):
    lumin_p = get_MC_luminosity(mcGen_files, period=i, charge="muPlus")
    lumin_m = get_MC_luminosity(mcGen_files, period=i, charge="muMinus") 
    print("Period:", i, ", MC_Lumin_p:", lumin_p)
    print("Period:", i, ", MC_Lumin_m:", lumin_m)
    intLumMC_muPlus.append(lumin_p)
    intLumMC_muMinus.append(lumin_m)

  # Get the normalization constants
  const_muPlus = [l_data / l_mc for l_data, l_mc in zip(intLum_muPlus, intLumMC_muPlus)]
  const_muMinus = [l_data / l_mc for l_data, l_mc in zip(intLum_muMinus, intLumMC_muMinus)]
  print("Period:", i, ", Const_p:", const_muPlus[0])
  print("Period:", i, ", Const_m:", const_muMinus[0])

  # Plot the data for the exclusivity variable 
  plot_variable_allPeriods(const_muPlus, const_muMinus, scale_factor=5.00, n_bins=27) # default plot is nu (full range)
  plot_variable_allPeriods(const_muPlus, const_muMinus, var="Q2_fit", scale_factor=0.56, n_bins=22, x_min=1, x_max=11, logy=True, logx=True)

  # *** Save values *** 
  write_constants_BH(intLum_muPlus, intLum_muMinus, const_muPlus, const_muMinus)

if __name__ == "__main__":

  main()