import os
import glob
import ROOT

# **********************************
# Real data 
real_dir = "/Users/gursimran/cern/2016_data/real/"
real_files = (os.path.join(real_dir, "filtered_P04.root"),
              os.path.join(real_dir, "filtered_P05.root"),
              os.path.join(real_dir, "filtered_P06.root"),
              os.path.join(real_dir, "filtered_P07.root"),
              os.path.join(real_dir, "filtered_P08.root"),
              os.path.join(real_dir, "filtered_P09.root"))

tree_name = "USR970_filtered"
chain = ROOT.TChain(tree_name)
for f in real_files:
  chain.Add(f)

n_total = chain.GetEntries()
print(f"Processing {n_total} entries in tree '{tree_name}' from {len(real_files)} files...")


# **********************************
# Save histograms
def saveHistogramsToPdf(histograms, output_pdf):
  c0 = ROOT.TCanvas("c0", "pull distributions", 1200, 900)
  c0.Print(output_pdf + "[")

  for i in range(0, len(histograms), 4):
    c0.Clear()
    c0.Divide(2, 2)

    for pad, hist in enumerate(histograms[i:i + 4], start=1):
      c0.cd(pad)
      hist.Draw()

    c0.Update()
    c0.Print(output_pdf)

  c0.Print(output_pdf + "]")


# **********************************
# Make the pull distributions
def makePullDist():
  # Define histograms
  h_inMu_x = ROOT.TH1F("hReal_inMu_x", ";#DeltaX_{#mu} / #sigma_{X_{#mu}}", 25, -10, 10)
  h_inMu_y = ROOT.TH1F("hReal_inMu_y", ";#DeltaY_{#mu} / #sigma_{Y_{#mu}}", 25, -10, 10)
  h_inMu_px = ROOT.TH1F("hReal_inMu_px", ";#Delta(P_{#mu})_{X} / #sigma_{(P_{#mu})_{X}}", 25, -10, 10)
  h_inMu_py = ROOT.TH1F("hReal_inMu_py", ";#Delta(P_{#mu})_{Y} / #sigma_{(P_{#mu})_{Y}}", 25, -10, 10)
  h_inMu_pz = ROOT.TH1F("hReal_inMu_pz", ";#Delta(P_{#mu})_{Z} / #sigma_{(P_{#mu})_{Z}}", 25, -10, 10)

  h_outMu_x = ROOT.TH1F("hReal_outMu_x", ";#DeltaX_{#mu'} / #sigma_{X_{#mu'}}", 25, -10, 10)
  h_outMu_y = ROOT.TH1F("hReal_outMu_y", ";#DeltaY_{#mu'} / #sigma_{Y_{#mu'}}", 25, -10, 10)
  h_outMu_px = ROOT.TH1F("hReal_outMu_px", ";#Delta(P_{#mu'})_{X} / #sigma_{(P_{#mu'})_{X}}", 25, -10, 10)
  h_outMu_py = ROOT.TH1F("hReal_outMu_py", ";#Delta(P_{#mu'})_{Y} / #sigma_{(P_{#mu'})_{Y}}", 25, -10, 10)
  h_outMu_pz = ROOT.TH1F("hReal_outMu_pz", ";#Delta(P_{#mu'})_{Z} / #sigma_{(P_{#mu'})_{Z}}", 25, -10, 10)

  h_gamma_x = ROOT.TH1F("hReal_gamma_x", ";#DeltaX_{#gamma} / #sigma_{X_{#gamma}}", 25, -10, 10)
  h_gamma_y = ROOT.TH1F("hReal_gamma_y", ";#DeltaY_{#gamma} / #sigma_{Y_{#gamma}}", 25, -10, 10)
  h_gamma_E = ROOT.TH1F("hReal_gamma_E", ";#DeltaE_{#gamma} / #sigma_{E_{#gamma}}", 25, -10, 10)

  h_pro_p = ROOT.TH1F("hReal_pro_x", ";#Delta{P_{p}} / #sigma_{P_{p}}", 25, -10, 10)
  h_pro_zA = ROOT.TH1F("hReal_pro_zA", ";#DeltaZ_{A} / #sigma_{Z_{A}}", 25, -10, 10)
  h_pro_rA = ROOT.TH1F("hReal_pro_rA", ";#Deltar_{A} / #sigma_{r_{A}}", 25, -10, 10)
  h_pro_phiA = ROOT.TH1F("hReal_pro_phiA", ";#Delta#Phi_{A} / #sigma_{#Phi_{A}}", 25, -10, 10)
  h_pro_zB = ROOT.TH1F("hReal_pro_zB", ";#DeltaZ_{B} / #sigma_{Z_{B}}", 25, -10, 10)
  h_pro_rB = ROOT.TH1F("hReal_pro_rB", ";#Deltar_{B} / #sigma_{r_{B}}", 25, -10, 10)
  h_pro_phiB = ROOT.TH1F("hReal_pro_phiB", ";#Delta#Phi_{B} / #sigma_{#Phi_{B}}", 25, -10, 10)

  # Loop over events and fill histograms
  for event in chain: 
    nu = event.nu # using the unfitted version here, honestly idk what to use here 
    if 10 < nu < 32: 
      h_inMu_x.Fill((event.pVtx_vec.X() - event.pVtxFit_vec.X()) / event.inMu_sigmaX)
      h_inMu_y.Fill((event.pVtx_vec.Y() - event.pVtxFit_vec.Y()) / event.inMu_sigmaY)
      h_inMu_px.Fill((event.inMu_TL.Px() - event.inMuFit_TL.Px()) / event.outMu_sigmaPx)
      h_inMu_py.Fill((event.inMu_TL.Py() - event.inMuFit_TL.Py()) / event.inMu_sigmaPy)
      h_inMu_pz.Fill((event.inMu_TL.Pz() - event.inMuFit_TL.Pz()) / event.inMu_sigmaPz)

      h_outMu_x.Fill((event.outMu_vec.X() - event.outMuFit_vec.X()) / event.outMu_sigmaX)
      h_outMu_y.Fill((event.outMu_vec.Y() - event.outMuFit_vec.Y()) / event.outMu_sigmaY)
      h_outMu_px.Fill((event.outMu_TL.Px() - event.outMuFit_TL.Px()) / event.outMu_sigmaPx)
      h_outMu_py.Fill((event.outMu_TL.Py() - event.outMuFit_TL.Py()) / event.outMu_sigmaPy)
      h_outMu_pz.Fill((event.outMu_TL.Pz() - event.outMuFit_TL.Pz()) / event.outMu_sigmaPz)

      # cluster_TL = (x, y, z, E), clusterFit_vec = (x, y, E)
      h_gamma_x.Fill((event.cluster_TL.X() - event.clusterFit_vec.X()) / event.gamma_sigmaX)
      h_gamma_y.Fill((event.cluster_TL.Y() - event.clusterFit_vec.Y()) / event.gamma_sigmaY)
      h_gamma_E.Fill((event.cluster_TL.T() - event.clusterFit_vec.Z()) / event.gamma_sigmaE)

      h_pro_p.Fill((event.p_camera_TL.P() - event.protonFit_TL.P()) / event.proton_sigmaP)
      if event.ringA_sigmaZ != 0: # protect against division by zero 
        h_pro_zA.Fill((event.posRingA_vec.Z() - event.posRingAFit_vec.Z()) / event.ringA_sigmaZ)
      # r = sqrt(x^2 + y^2) -> radial vector in the XY plane 
      h_pro_rA.Fill((event.posRingA_vec.Perp() - event.posRingAFit_vec.Perp()) / event.ringA_sigmaR)
      h_pro_phiA.Fill((event.posRingA_vec.Phi() - event.posRingAFit_vec.Phi()) / event.ringA_sigmaPhi)
      if event.ringB_sigmaZ != 0: # protect against division by zero 
        h_pro_zB.Fill((event.posRingB_vec.Z() - event.posRingBFit_vec.Z()) / event.ringB_sigmaZ)
      # r = sqrt(x^2 + y^2) -> radial vector in the XY plane 
      h_pro_rB.Fill((event.posRingB_vec.Perp() - event.posRingBFit_vec.Perp()) / event.ringB_sigmaR)
      h_pro_phiB.Fill((event.posRingB_vec.Phi() - event.posRingBFit_vec.Phi()) / event.ringB_sigmaPhi)

  # Save the histograms to a multi-page PDF, four plots per page
  inMu_hists = (h_inMu_x, h_inMu_y, h_inMu_px, h_inMu_py, h_inMu_pz)
  outMu_hists = (h_outMu_x, h_outMu_y, h_outMu_px, h_outMu_py, h_outMu_pz,)
  gamma_hists = (h_gamma_x, h_gamma_y, h_gamma_E,)
  pro_hists = (h_pro_p, h_pro_zA, h_pro_rA, h_pro_phiA, h_pro_zB, h_pro_rB, h_pro_phiB)

  saveHistogramsToPdf(inMu_hists, "inMu_pulls.pdf")
  #saveHistogramsToPdf(outMu_hists, "outMu_pulls.pdf")
  #saveHistogramsToPdf(gamma_hists, "gamma_pulls.pdf")
  #saveHistogramsToPdf(pro_hists, "proton_pulls.pdf")


# **********************************
# Main function 
def main():
  makePullDist()

if __name__ == "__main__":
  main()
