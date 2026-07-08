import os
import math
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# *******************************************************************
# *                    *** DATA PREP ***                            *
# *******************************************************************
# **********************************
# Output directory 
out_dir = "/Users/gursimran/cern/kinFitPlots"

# **********************************
# Real data 
real_dir = "/Users/gursimran/cern/2016_data/real/"
real_files = (os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root"),)

chain_real = ROOT.TChain("USR970_filtered")
for f in real_files:
  chain_real.Add(f)

n_total = chain_real.GetEntries()
print(f"Processing {n_total} real entries in tree 'USR970_filtered' from {len(real_files)} files...")

# **********************************
# BH MC Data
hepBH_dir = "/Users/gursimran/cern/2016_data/BH"
hepBH_files = (os.path.join(hepBH_dir, "filtered_P04_muPlus.root"), os.path.join(hepBH_dir, "filtered_P04_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P05_muPlus.root"), os.path.join(hepBH_dir, "filtered_P05_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P06_muPlus.root"), os.path.join(hepBH_dir, "filtered_P06_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P07_muPlus.root"), os.path.join(hepBH_dir, "filtered_P07_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P08_muPlus.root"), os.path.join(hepBH_dir, "filtered_P08_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P09_muPlus.root"), os.path.join(hepBH_dir, "filtered_P09_muMinus.root"),)

chain_hepBH = ROOT.TChain("USR970_filtered")
for f in hepBH_files:
  chain_hepBH.Add(f)

n_total = chain_hepBH.GetEntries()
print(f"Processing {n_total} Hepgen_BH entries in tree 'USR970_filtered' from {len(hepBH_files)} files...")


# **********************************
# Hepgen pi0 Data 
hepPi0_dir = "/Users/gursimran/cern/2016_data/HepgenPi0"
hepPi0_files = (os.path.join(hepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(hepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(hepPi0_dir, "filtered_P09_muMinus.root"),)

chain_hepPi0 = ROOT.TChain("USR970_filtered")
for f in hepPi0_files:
  chain_hepPi0.Add(f)

n_total = chain_hepPi0.GetEntries()
print(f"Processing {n_total} Hepgen_Pi0 entries in tree 'USR970_filtered' from {len(hepPi0_files)} files...")


# **********************************
# Lepto pi0 Data 
lepPi0_dir = "/Users/gursimran/cern/2016_data/LeptoPi0"
lepPi0_files = (os.path.join(lepPi0_dir, "filtered_P04_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P04_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P05_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P05_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P06_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P06_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P07_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P07_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P08_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P08_muMinus.root"),
                os.path.join(lepPi0_dir, "filtered_P09_muPlus.root"), os.path.join(lepPi0_dir, "filtered_P09_muMinus.root"),)

chain_lepPi0 = ROOT.TChain("USR970_filtered")
for f in lepPi0_files:
  chain_lepPi0.Add(f)

n_total = chain_lepPi0.GetEntries()
print(f"Processing {n_total} Lepto_Pi0 entries in tree 'USR970_filtered' from {len(lepPi0_files)} files...")


# *******************************************************************
# *                 *** PLOT CONFIGS AND HELPERS ***                *
# *******************************************************************
# **********************************
# Extra params in "conf_level" tuples are for xMin and xMax - plots have variable ranges unlike pulls 
PLOT_CONFIGS = {
  "inMu": (
    ("inMu_x", "#DeltaX_{#mu} / #sigma_{X_{#mu}}"),
    ("inMu_y", "#DeltaY_{#mu} / #sigma_{Y_{#mu}}"),
    ("inMu_px", "#Delta(P_{#mu})_{X} / #sigma_{(P_{#mu})_{X}}"),
    ("inMu_py", "#Delta(P_{#mu})_{Y} / #sigma_{(P_{#mu})_{Y}}"),
    ("inMu_pz", "#Delta(P_{#mu})_{Z} / #sigma_{(P_{#mu})_{Z}}"),
  ),
  "outMu": (
    ("outMu_x", "#DeltaX_{#mu'} / #sigma_{X_{#mu'}}"),
    ("outMu_y", "#DeltaY_{#mu'} / #sigma_{Y_{#mu'}}"),
    ("outMu_px", "#Delta(P_{#mu'})_{X} / #sigma_{(P_{#mu'})_{X}}"),
    ("outMu_py", "#Delta(P_{#mu'})_{Y} / #sigma_{(P_{#mu'})_{Y}}"),
    ("outMu_pz", "#Delta(P_{#mu'})_{Z} / #sigma_{(P_{#mu'})_{Z}}"),
  ),
  "gamma": (
    ("gamma_x", "#DeltaX_{#gamma} / #sigma_{X_{#gamma}}"),
    ("gamma_y", "#DeltaY_{#gamma} / #sigma_{Y_{#gamma}}"),
    ("gamma_E", "#DeltaE_{#gamma} / #sigma_{E_{#gamma}}"),
  ),
  "proton": (
    ("pro_p", "#Delta{P_{p}} / #sigma_{P_{p}}"),
    ("ringA_z", "#DeltaZ_{A} / #sigma_{Z_{A}}"),
    ("ringA_r", "#Deltar_{A} / #sigma_{r_{A}}"),
    ("ringA_phi", "#Delta#Phi_{A} / #sigma_{#Phi_{A}}"),
    ("ringB_z", "#DeltaZ_{B} / #sigma_{Z_{B}}"),
    ("ringB_r", "#Deltar_{B} / #sigma_{r_{B}}"),
    ("ringB_phi", "#Delta#Phi_{B} / #sigma_{#Phi_{B}}"),
  ),
  "conf_level": (
    ("cl", "Confidence level", 0, 1),
    ("chi2", "#chi^{2}", 0, 70),
    ("red_chi2", "#chi^{2}/ndf", 0, 10),
    ("cons_mom", "p_{#mu} + p_{p} - p_{#mu'} - p_{#gamma} - p_{p'}", -1, 1),
    ("worst_pull", "Largest pull^{2} contributor", 0, 20),
  )
}

# **********************************
INMU_PULL_EXPRESSIONS = {
  "inMu_x": "(pVtx_vec.X() - pVtxFit_vec.X()) / inMu_sigmaX",
  "inMu_y": "(pVtx_vec.Y() - pVtxFit_vec.Y()) / inMu_sigmaY",
  "inMu_px": "(inMu_TL.Px() - inMuFit_TL.Px()) / inMu_sigmaPx",
  "inMu_py": "(inMu_TL.Py() - inMuFit_TL.Py()) / inMu_sigmaPy",
  "inMu_pz": "(inMu_TL.Pz() - inMuFit_TL.Pz()) / inMu_sigmaPz",
}

OUTMU_PULL_EXPRESSIONS = {
  "outMu_x": "(outMu_vec.X() - outMuFit_vec.X()) / outMu_sigmaX",
  "outMu_y": "(outMu_vec.Y() - outMuFit_vec.Y()) / outMu_sigmaY",
  "outMu_px": "(outMu_TL.Px() - outMuFit_TL.Px()) / outMu_sigmaPx",
  "outMu_py": "(outMu_TL.Py() - outMuFit_TL.Py()) / outMu_sigmaPy",
  "outMu_pz": "(outMu_TL.Pz() - outMuFit_TL.Pz()) / outMu_sigmaPz",
}

GAMMA_PULL_EXPRESSIONS = {
  "gamma_x": "(cluster_TL.X() - clusterFit_vec.X()) / gamma_sigmaX",
  "gamma_y": "(cluster_TL.Y() - clusterFit_vec.Y()) / gamma_sigmaY",
  "gamma_E": "(gamma_TL.E() - gammaFit_TL.E()) / gamma_sigmaE",
}

PROTON_PULL_EXPRESSIONS = {
  "pro_p": "(p_camera_TL.P() - protonFit_TL.P()) / proton_sigmaP",
  "ringA_z": "(posRingA_vec.Z() - posRingAFit_vec.Z()) / ringA_sigmaZ",
  "ringA_r": "(posRingA_vec.Perp() - posRingAFit_vec.Perp()) / ringA_sigmaR",
  "ringA_phi": "(posRingA_vec.Phi() - posRingAFit_vec.Phi()) / ringA_sigmaPhi",
  "ringB_z": "(posRingB_vec.Z() - posRingBFit_vec.Z()) / ringB_sigmaZ",
  "ringB_r": "(posRingB_vec.Perp() - posRingBFit_vec.Perp()) / ringB_sigmaR",
  "ringB_phi": "(posRingB_vec.Phi() - posRingBFit_vec.Phi()) / ringB_sigmaPhi",
}

PULL_EXPRESSIONS = {
  "inMu": INMU_PULL_EXPRESSIONS,
  "outMu": OUTMU_PULL_EXPRESSIONS,
  "gamma": GAMMA_PULL_EXPRESSIONS,
  "proton": PROTON_PULL_EXPRESSIONS,
}

PULL_NAMES = [
  name
  for expressions in PULL_EXPRESSIONS.values()
  for name in expressions
]

# **********************************
CL_EXPRESSIONS = {
  "cl": "TMath::Prob(chi2_fit, ndf_fit)",
  "chi2": "chi2_fit",
  "red_chi2": "chi2_fit / ndf_fit",
  "cons_mom": "inMuFit_TL.P() + targetFit_TL.P() - outMuFit_TL.P() - gammaFit_TL.P() - protonFit_TL.P()", 
}

# **********************************
DATA_CONFIGS = {
  "real": {
    "chain": chain_real,
    "dataType": "real",
    "label": "Real",
    "lineColor": ROOT.kBlack,
    "fillColor": ROOT.kGray+1,
    "markerStyle": 33,
  },
  "hep_BH": {
    "chain": chain_hepBH,
    "dataType": "Hepgen_BH",
    "label": "Hepgen BH",
    "lineColor": ROOT.kRed,
    "fillColor": ROOT.kRed,
    "markerStyle": None,
  },
  "hep_pi0": {
    "chain": chain_hepPi0,
    "dataType": "Hepgen_Pi0",
    "label": "Hepgen Pi0",
    "lineColor": ROOT.kBlue,
    "fillColor": ROOT.kBlue,
    "markerStyle": None,
  },
  "lep_pi0": {
    "chain": chain_lepPi0,
    "dataType": "Lepto_Pi0",
    "label": "Lepto Pi0",
    "lineColor": ROOT.kAzure + 7,
    "fillColor": ROOT.kAzure + 7,
    "markerStyle": None,
  },
}

# **********************************
# Normalize a histogram so its total area is 1 for plotting 
def normalizeHist(hist):
  integral = hist.Integral()
  if integral > 0:
    hist.Scale(1.0 / integral)

# **********************************
# Reformat the x and y axes for the plot 
def formatAxes(hist):
  x_axis = hist.GetXaxis()
  y_axis = hist.GetYaxis()

  x_axis.SetTitleSize(0.04)
  y_axis.SetTitleSize(0.04)
  x_axis.SetLabelSize(0.035)
  y_axis.SetLabelSize(0.035)
  x_axis.SetTitleOffset(1.25)
  y_axis.SetTitleOffset(1.25)

# **********************************
# Choose which datasets to process
def selectDatasets(data_set="all"):
  if data_set == "all":
    return tuple(DATA_CONFIGS)

  if isinstance(data_set, str):
    selected = (data_set,)
  else:
    selected = tuple(data_set)

  if len(selected) == 0:
    raise ValueError(f"No datasets selected. Options are: all, {', '.join(DATA_CONFIGS)}")

  unknown = [dataset for dataset in selected if dataset not in DATA_CONFIGS]
  if unknown:
    raise ValueError(f"Unknown dataset(s): {', '.join(unknown)}. Options are: all, {', '.join(DATA_CONFIGS)}")

  return selected

# **********************************
# Build selection string for TChain::Draw
def buildSelection(selection="nu > 10 && nu < 32", clCut=None):
  if clCut is None:
    return selection

  if clCut < 0 or clCut > 1:
    raise ValueError("clCut must be between 0 and 1")

  return f"({selection}) && (TMath::Prob(chi2_fit, ndf_fit) > {clCut})"

# **********************************
# Label for CL-cut output directory
def formatClCutLabel(clCut):
  return str(clCut).replace(".", "p")


# *******************************************************************
# *                   *** PULL HELPERS ***                          *
# *******************************************************************
# **********************************
# Initialize the histograms  
def initPullHists(dataType="real", groupName="inMu"):
  if groupName == "all":
    groups = (
      group
      for name, group in PLOT_CONFIGS.items()
      if name != "conf_level"
    )
  else:
    if groupName not in PLOT_CONFIGS:
      raise ValueError(f"Unknown histogram group '{groupName}'. Options are: all, {', '.join(PLOT_CONFIGS)}")
    groups = (PLOT_CONFIGS[groupName],)

  histograms = {}
  for group in groups:
    for name, title in group:
      histograms[name] = ROOT.TH1F(f"h{dataType}_{name}", f";{title}; Normalized events", 35, -10, 10)
  return histograms

# **********************************
# Fill the histograms using Draw function
def fillPullHistsWithDraw(chain, histograms, groupName="inMu", selection="nu > 10 && nu < 32"):
  if groupName == "all":
    expression_groups = PULL_EXPRESSIONS.values()
  else:
    if groupName not in PULL_EXPRESSIONS:
      raise ValueError(f"Unknown pull group '{groupName}'. Options are: all, {', '.join(PULL_EXPRESSIONS)}")
    expression_groups = (PULL_EXPRESSIONS[groupName],)

  for expressions in expression_groups:
    for name, expression in expressions.items():
      if name not in histograms:
        continue

      hist = histograms[name]
      hist.Reset()
      chain.Draw(f"{expression} >> {hist.GetName()}", selection, "goff")

# **********************************
# Format and save each plot as a png 
def savePullHistsToPngs(histsByDataset, outputDir=out_dir):
  os.makedirs(outputDir, exist_ok=True)

  c0 = ROOT.TCanvas("c0", "histogram canvas", 1000, 700)
  c0.SetLeftMargin(0.12)
  c0.SetRightMargin(0.05)
  c0.SetBottomMargin(0.14)
  c0.SetTopMargin(0.07)

  normalize_plots = len(histsByDataset) > 1
  first_dataset = next(iter(histsByDataset))
  for name in histsByDataset[first_dataset]:
    if any(name not in hists for hists in histsByDataset.values()):
      continue

    if all(hists[name].GetEntries() == 0 for hists in histsByDataset.values()):
      continue

    plot_hists = {}
    for dataset, hists in histsByDataset.items():
      hist_plot = hists[name].Clone(f"{hists[name].GetName()}_plot")
      hist_plot.SetDirectory(0)
      if normalize_plots:
        normalizeHist(hist_plot)

      style = DATA_CONFIGS[dataset]
      hist_plot.SetLineColor(style["lineColor"])
      hist_plot.SetLineWidth(2)
      if style["fillColor"] is not None:
        hist_plot.SetFillColorAlpha(style["fillColor"], 0.2)
      if style["markerStyle"] is not None:
        hist_plot.SetMarkerColor(style["lineColor"])
        hist_plot.SetMarkerStyle(style["markerStyle"])
        hist_plot.SetMarkerSize(2.0)

      plot_hists[dataset] = hist_plot

    y_max = max(hist.GetMaximum() for hist in plot_hists.values())
    axis_dataset = next((dataset for dataset in plot_hists if dataset != "real"), next(iter(plot_hists)))
    axis_hist = plot_hists[axis_dataset]
    if y_max > 0:
      axis_hist.SetMaximum(1.2 * y_max)
    axis_hist.SetMinimum(0)
    axis_hist.GetYaxis().SetTitle("Normalized events" if normalize_plots else "Events")
    formatAxes(axis_hist)

    fit_real = None
    fit_text = None
    if "real" in plot_hists and plot_hists["real"].GetEntries() > 0:
      hist_real_plot = plot_hists["real"]
      fit_real = ROOT.TF1(f"fit_{name}", "gaus", -4, 4)
      fit_real.SetLineColor(ROOT.kBlack)
      fit_real.SetLineStyle(2)
      fit_real.SetLineWidth(2)
      hist_real_plot.Fit(fit_real, "RQ0")

      fit_text = ROOT.TPaveText(0.16, 0.68, 0.40, 0.86, "NDC")
      fit_text.SetFillColor(0)
      fit_text.SetFillStyle(1001)
      fit_text.SetBorderSize(1)
      fit_text.SetTextAlign(12)
      fit_text.SetTextSize(0.035)
      fit_text.AddText("Real Gaussian Fit")
      fit_text.AddText(f"#mu = {fit_real.GetParameter(1):.3f}")
      fit_text.AddText(f"#sigma = {fit_real.GetParameter(2):.3f}")
      if fit_real.GetNDF() > 0:
        fit_text.AddText(f"#chi^{{2}}/ndf = {fit_real.GetChisquare() / fit_real.GetNDF():.2f}")

    c0.Clear()
    axis_hist.Draw("E" if axis_dataset == "real" else "HIST")
    for dataset, hist_plot in plot_hists.items():
      if dataset == axis_dataset:
        continue

      draw_option = "E SAME" if dataset == "real" else "HIST SAME"
      hist_plot.Draw(draw_option)

    if fit_real is not None:
      fit_real.Draw("SAME")
    if fit_text is not None:
      fit_text.Draw()

    legend = ROOT.TLegend(0.60, 0.70, 0.90, 0.90)
    legend.SetTextSize(0.035)
    if "real" in plot_hists:
      legend.AddEntry(plot_hists["real"], DATA_CONFIGS["real"]["label"], "lep")
    if fit_real is not None:
      legend.AddEntry(fit_real, "Real Gaussian fit", "l")
    for dataset, hist_plot in plot_hists.items():
      if dataset == "real":
        continue

      legend.AddEntry(hist_plot, DATA_CONFIGS[dataset]["label"], "lf")
    legend.Draw()

    c0.Update()
    c0.SaveAs(os.path.join(outputDir, f"{name}.png"))


# *******************************************************************
# *              *** CONFIDENCE LEVEL HELPERS ***                   *
# *******************************************************************
# **********************************
# Initialize confidence level histograms
def initConfHists(dataType="real"):
  histograms = {}
  for name, title, x_min, x_max in PLOT_CONFIGS["conf_level"]:
    if name == "worst_pull":
      histograms[name] = ROOT.TH1F(f"h{dataType}_{name}", f";{title};Normalized events", len(PULL_NAMES), x_min, x_max)
      for i, pull_name in enumerate(PULL_NAMES, start=1):
        histograms[name].GetXaxis().SetBinLabel(i, pull_name)
    else:
      histograms[name] = ROOT.TH1F(f"h{dataType}_{name}", f";{title};Normalized events", 35, x_min, x_max)
  return histograms

# **********************************
# Fill confidence-level histograms
def fillConfHistsWithDraw(chain, histograms, selection="nu > 10 && nu < 32"):
  for name, expression in CL_EXPRESSIONS.items():
    if name not in histograms:
      continue

    hist = histograms[name]
    hist.Reset()
    chain.Draw(f"{expression} >> {hist.GetName()}", selection, "goff")

  if "worst_pull" in histograms:
    fillWorstPullHist(chain, histograms["worst_pull"], selection)

# **********************************
# Fill the largest pull contribution to chi2
def fillWorstPullHist(chain, hist, selection="nu > 10 && nu < 32"):
  hist.Reset()

  selection_formula = ROOT.TTreeFormula("worst_pull_selection", selection, chain)
  pull_formulas = {
    name: ROOT.TTreeFormula(f"worst_pull_{name}", expression, chain)
    for name, expression in (
      (pull_name, pull_expression)
      for expressions in PULL_EXPRESSIONS.values()
      for pull_name, pull_expression in expressions.items()
    )
  }

  current_tree_number = -1
  for i_entry in range(chain.GetEntries()):
    chain.GetEntry(i_entry)
    if chain.GetTreeNumber() != current_tree_number:
      current_tree_number = chain.GetTreeNumber()
      selection_formula.UpdateFormulaLeaves()
      for formula in pull_formulas.values():
        formula.UpdateFormulaLeaves()

    if not selection_formula.EvalInstance():
      continue

    contribs = {}
    for name, formula in pull_formulas.items():
      pull = formula.EvalInstance()
      if not math.isfinite(pull):
        continue
      contribs[name] = pull * pull

    if len(contribs) == 0:
      continue

    worst_name = max(contribs, key=contribs.get)
    hist.Fill(PULL_NAMES.index(worst_name))

# **********************************
# Save confidence-level plots
def saveConfLevelToPngs(histsByDataset, outputDir=out_dir):
  os.makedirs(outputDir, exist_ok=True)

  c0 = ROOT.TCanvas("c_conf", "confidence canvas", 1000, 700)
  c0.SetLeftMargin(0.12)
  c0.SetRightMargin(0.05)
  c0.SetBottomMargin(0.14)
  c0.SetTopMargin(0.07)

  normalize_plots = len(histsByDataset) > 1
  first_dataset = next(iter(histsByDataset))
  for name in histsByDataset[first_dataset]:
    if any(name not in hists for hists in histsByDataset.values()):
      continue

    if all(hists[name].GetEntries() == 0 for hists in histsByDataset.values()):
      continue

    plot_hists = {}
    for dataset, hists in histsByDataset.items():
      hist_plot = hists[name].Clone(f"{hists[name].GetName()}_plot")
      hist_plot.SetDirectory(0)
      if normalize_plots:
        normalizeHist(hist_plot)

      style = DATA_CONFIGS[dataset]
      hist_plot.SetLineColor(style["lineColor"])
      hist_plot.SetLineWidth(2)
      if style["fillColor"] is not None:
        hist_plot.SetFillColorAlpha(style["fillColor"], 0.2)
      if style["markerStyle"] is not None:
        hist_plot.SetMarkerColor(style["lineColor"])
        hist_plot.SetMarkerStyle(style["markerStyle"])
        hist_plot.SetMarkerSize(2.0)

      plot_hists[dataset] = hist_plot

    y_max = max(hist.GetMaximum() for hist in plot_hists.values())
    axis_dataset = next((dataset for dataset in plot_hists if dataset != "real"), next(iter(plot_hists)))
    axis_hist = plot_hists[axis_dataset]
    if y_max > 0:
      axis_hist.SetMaximum(1.2 * y_max)
    axis_hist.SetMinimum(0)
    axis_hist.GetYaxis().SetTitle("Normalized events" if normalize_plots else "Events")
    formatAxes(axis_hist)

    c0.Clear()
    axis_hist.Draw("HIST")
    for dataset, hist_plot in plot_hists.items():
      if dataset == axis_dataset:
        continue

      draw_option = "HIST SAME"
      hist_plot.Draw(draw_option)

    legend = ROOT.TLegend(0.72, 0.76, 0.92, 0.90)
    legend.SetTextSize(0.03)
    for dataset, hist_plot in plot_hists.items():
      legend_option = "l" if DATA_CONFIGS[dataset]["fillColor"] is None else "lf"
      legend.AddEntry(hist_plot, DATA_CONFIGS[dataset]["label"], legend_option)
    legend.Draw()

    c0.Update()
    c0.SaveAs(os.path.join(outputDir, f"conf_{name}.png"))


# *******************************************************************
# *                    *** WRAPPER FUNCTIONS ***                    *
# *******************************************************************
# **********************************
# Wrapper to make the pull distributions
# Can specify what data_set to use and which set of pulls to make 
# ex. makePulls(data_set=["hep_pi0", "lep_pi0"], groupName="inMu")
# ex. makePulls(data_set="real", groupName="gamma", clCut=0.1)
def makePulls(data_set="all", groupName="all", selection="nu > 10 && nu < 32", clCut=None, outputDir=None):
  if groupName != "all" and groupName not in PULL_EXPRESSIONS:
    raise ValueError(f"Unknown pull group '{groupName}'. Options are: all, {', '.join(PULL_EXPRESSIONS)}")

  selection = buildSelection(selection=selection, clCut=clCut)
  if outputDir is None:
    outputDir = out_dir
    if clCut is not None:
      outputDir = os.path.join(out_dir, f"CLgt{formatClCutLabel(clCut)}")

  selected_datasets = selectDatasets(data_set)
  hists_by_dataset = {}
  for dataset in selected_datasets:
    dataset_config = DATA_CONFIGS[dataset]
    histograms = initPullHists(dataType=dataset_config["dataType"], groupName=groupName)
    fillPullHistsWithDraw(dataset_config["chain"], histograms, groupName=groupName, selection=selection)
    hists_by_dataset[dataset] = histograms

  savePullHistsToPngs(hists_by_dataset, outputDir=outputDir)

# **********************************
# Wrapper to make the confidence level distributions 
def makeConfLevel(data_set="real"):
  selected_datasets = selectDatasets(data_set)
  hists_by_dataset = {}
  for dataset in selected_datasets:
    dataset_config = DATA_CONFIGS[dataset]
    histograms = initConfHists(dataType=dataset_config["dataType"])
    fillConfHistsWithDraw(dataset_config["chain"], histograms)
    hists_by_dataset[dataset] = histograms

  saveConfLevelToPngs(hists_by_dataset)


# *******************************************************************
# *                       *** MAIN ***                              *
# *******************************************************************
# **********************************
# Main function - exclude or include functions here 
def main():
  #makePulls(data_set="all", groupName="all")
  #makeConfLevel(data_set="real")
  makePulls(data_set="all", groupName="gamma", clCut=0.1)

if __name__ == "__main__":
  main()
