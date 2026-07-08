import os
import ROOT

ROOT.gROOT.SetBatch(True)

# **********************************
# Real data 
real_dir = "/Users/gursimran/cern/2016_data/real/"
real_files = (os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root"))

chain = ROOT.TChain("USR970_filtered")
for f in real_files:
  chain.Add(f)

n_total = chain.GetEntries()
print(f"Processing {n_total} entries in tree 'USR970_filtered' from {len(real_files)} files...")

# **********************************
# BH MC Data
hepBH_dir = "/Users/gursimran/cern/2016_data/BH"
""" hepBH_files = (os.path.join(hepBH_dir, "filtered_P04_muPlus.root"), os.path.join(hepBH_dir, "filtered_P04_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P05_muPlus.root"), os.path.join(hepBH_dir, "filtered_P05_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P06_muPlus.root"), os.path.join(hepBH_dir, "filtered_P06_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P07_muPlus.root"), os.path.join(hepBH_dir, "filtered_P07_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P08_muPlus.root"), os.path.join(hepBH_dir, "filtered_P08_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P09_muPlus.root"), os.path.join(hepBH_dir, "filtered_P09_muMinus.root"),) """

""" hepBH_files = (os.path.join(hepBH_dir, "filtered_P04_muPlus.root"), os.path.join(hepBH_dir, "filtered_P04_muMinus.root"),
               os.path.join(hepBH_dir, "filtered_P05_muPlus.root"), os.path.join(hepBH_dir, "filtered_P05_muMinus.root"),)

chain_hepBH = ROOT.TChain("USR970_filtered")
for f in hepBH_files:
  chain_hepBH.Add(f)

n_total = chain_hepBH.GetEntries()
print(f"Processing {n_total} entries in tree 'USR970_filtered' from {len(hepBH_files)} files...") """


# **********************************
# Save histograms
def saveHistogramsToPngs(histograms, output_dir="kinFitPlots"):
  os.makedirs(output_dir, exist_ok=True)

  c0 = ROOT.TCanvas("c0", "histogram canvas", 800, 600)
  for hist in histograms:
    c0.Clear()
    hist.Draw()
    c0.Update()

    output_name = hist.GetName()
    if output_name.startswith("hReal_"):
      output_name = output_name[len("hReal_"):]
    elif output_name.startswith("h_"):
      output_name = output_name[len("h_"):]

    c0.SaveAs(os.path.join(output_dir, f"{output_name}.png"))


# **********************************
# Calculate single pull 
def getPull(meas, fit, sigma):
  if sigma == 0:
    return None
  return (meas - fit) / sigma

# Calculate all pulls for an event
def getPulls(event):
  return {
    "inMu_x": getPull(event.pVtx_vec.X(), event.pVtxFit_vec.X(), event.inMu_sigmaX),
    "inMu_y": getPull(event.pVtx_vec.Y(), event.pVtxFit_vec.Y(), event.inMu_sigmaY),
    "inMu_px": getPull(event.inMu_TL.Px(), event.inMuFit_TL.Px(), event.inMu_sigmaPx),
    "inMu_py": getPull(event.inMu_TL.Py(), event.inMuFit_TL.Py(), event.inMu_sigmaPy),
    "inMu_pz": getPull(event.inMu_TL.Pz(), event.inMuFit_TL.Pz(), event.inMu_sigmaPz),

    "outMu_x": getPull(event.outMu_vec.X(), event.outMuFit_vec.X(), event.outMu_sigmaX),
    "outMu_y": getPull(event.outMu_vec.Y(), event.outMuFit_vec.Y(), event.outMu_sigmaY),
    "outMu_px": getPull(event.outMu_TL.Px(), event.outMuFit_TL.Px(), event.outMu_sigmaPx),
    "outMu_py": getPull(event.outMu_TL.Py(), event.outMuFit_TL.Py(), event.outMu_sigmaPy),
    "outMu_pz": getPull(event.outMu_TL.Pz(), event.outMuFit_TL.Pz(), event.outMu_sigmaPz),

    "gamma_x": getPull(event.cluster_TL.X(), event.clusterFit_vec.X(), event.gamma_sigmaX),
    "gamma_y": getPull(event.cluster_TL.Y(), event.clusterFit_vec.Y(), event.gamma_sigmaY),
    "gamma_E": getPull(event.gamma_TL.E(), event.gammaFit_TL.E(), event.gamma_sigmaE),

    "pro_p": getPull(event.p_camera_TL.P(), event.protonFit_TL.P(), event.proton_sigmaP),
    "ringA_z": getPull(event.posRingA_vec.Z(), event.posRingAFit_vec.Z(), event.ringA_sigmaZ),
    "ringA_r": getPull(event.posRingA_vec.Perp(), event.posRingAFit_vec.Perp(), event.ringA_sigmaR),
    "ringA_phi": getPull(event.posRingA_vec.Phi(), event.posRingAFit_vec.Phi(), event.ringA_sigmaPhi),
    "ringB_z": getPull(event.posRingB_vec.Z(), event.posRingBFit_vec.Z(), event.ringB_sigmaZ),
    "ringB_r": getPull(event.posRingB_vec.Perp(), event.posRingBFit_vec.Perp(), event.ringB_sigmaR),
    "ringB_phi": getPull(event.posRingB_vec.Phi(), event.posRingBFit_vec.Phi(), event.ringB_sigmaPhi),
  }

# **********************************
# Make the pull distributions
PULL_HIST_GROUPS = {
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
}

PULL_NAMES = [
  name
  for group in PULL_HIST_GROUPS.values()
  for name, _ in group
]


def makePullDist(chain):
  histograms = {}
  for group in PULL_HIST_GROUPS.values():
    for name, title in group:
      histograms[name] = ROOT.TH1F(f"hReal_{name}", f";{title};Events", 25, -10, 10)

  for event in chain:
    pulls = getPulls(event)
    for name, pull in pulls.items():
      if pull is not None:
        histograms[name].Fill(pull)

  for group_name, group in PULL_HIST_GROUPS.items():
    group_hists = tuple(histograms[name] for name, _ in group)
    saveHistogramsToPngs(group_hists)


# **********************************
# Make the confidence level and related distributions
def confLevel():
  h_cl = ROOT.TH1F("h_fitCL", ";Confidence level;Events", 50, 0, 1)
  h_chi2 = ROOT.TH1F("h_chi2", ";#chi^{2};Events", 100, 0, 70)
  h_ndf = ROOT.TH1F("h_ndf", ";ndf;Events", 20, 0, 20)
  h_chi2ndf = ROOT.TH1F("h_chi2ndf", ";#chi^{2}/ndf;Events", 100, 0, 10)

  h_consMom = ROOT.TH1F("h_consMom", ";p_{#mu} + p_{p} - p_{#mu'} - p_{#gamma} - p_{p'}", 50, -1, 1)
  h_consMom_lowCL = ROOT.TH1F("h_consMom_lowCL", "CL < 0.1;p_{#mu} + p_{p} - p_{#mu'} - p_{#gamma} - p_{p'}",50, -1, 1)
  h_consMom_highCL = ROOT.TH1F("h_consMom_highCL", "CL >= 0.1;p_{#mu} + p_{p} - p_{#mu'} - p_{#gamma} - p_{p'}",50, -1, 1)

  h_worst = ROOT.TH1F("h_worst_lowCL", "CL < 0.1;Largest pull^{2} contributor", len(PULL_NAMES), 0, len(PULL_NAMES))

  # Define additional histograms for the worst offender 
  h_gammaE_pull_lowCL = ROOT.TH1F("h_gammaE_pull_lowCL", "CL < 0.1;#Delta E_{#gamma}/#sigma_{E};Events", 25, -10, 10)
  h_gammaE_sigma_lowCL = ROOT.TH1F("h_gammaE_sigma_lowCL", "CL < 0.1;#sigma_{E_{#gamma}};Events", 25, 0, 2)
  h_gammaE_delta_lowCL = ROOT.TH1F("h_gammaE_delta_lowCL", "CL < 0.1;E_{#gamma}^{meas} - E_{#gamma}^{fit};Events", 25, -10, 10)

  for i, name in enumerate(PULL_NAMES, start=1):
    h_worst.GetXaxis().SetBinLabel(i, name)

  for event in chain:
    if event.nu > 10 and event.nu < 32: # fit just in the extraction region 
      chi2 = event.chi2_fit
      ndf = int(event.ndf_fit)

      cl = ROOT.TMath.Prob(chi2, ndf)
      chi2ndf = chi2 / ndf

      h_cl.Fill(cl)
      h_chi2.Fill(chi2)
      h_ndf.Fill(ndf)
      h_chi2ndf.Fill(chi2ndf)

      check_consMom = (event.inMuFit_TL.P() + event.targetFit_TL.P() - 
                       event.outMuFit_TL.P() - event.gammaFit_TL.P() - event.protonFit_TL.P())
      h_consMom.Fill(check_consMom)

      if cl < 0.1:
        h_consMom_lowCL.Fill(check_consMom)

        pulls = getPulls(event)
        contribs = {
          name: pull * pull
          for name, pull in pulls.items()
          if pull is not None
        }

        if len(contribs) == 0:
          continue

        worst_name = max(contribs, key=contribs.get)
        worst_pull = pulls[worst_name]

        h_worst.Fill(PULL_NAMES.index(worst_name))
        #print("CL =", cl, "chi2 =", chi2, "worst =", worst_name, "pull =", worst_pull)

        # Plot additional distribution for the worst offender
        deltaE = event.gamma_TL.E() - event.gammaFit_TL.Z()
        pullE = deltaE / event.gamma_sigmaE

        h_gammaE_pull_lowCL.Fill(pullE)
        h_gammaE_sigma_lowCL.Fill(event.gamma_sigmaE)
        h_gammaE_delta_lowCL.Fill(deltaE)

      else:
        h_consMom_highCL.Fill(check_consMom)

  other_hists = (h_cl, h_chi2, h_chi2ndf, h_ndf, h_consMom, h_consMom_lowCL, h_consMom_highCL, h_worst,
                 h_gammaE_pull_lowCL, h_gammaE_sigma_lowCL, h_gammaE_delta_lowCL)
  saveHistogramsToPngs(other_hists)



# **********************************
# Main function 
def main():
  makePullDist(chain_hepBH)
  #confLevel()

if __name__ == "__main__":
  main()
