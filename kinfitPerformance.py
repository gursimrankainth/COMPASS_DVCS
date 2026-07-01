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
  # outMu pulls
  h_outMu_x = ROOT.TH1F("hReal_outMu_x", "; #DeltaX_{#mu^{#prime}} / #sigma_{X_{#mu^{#prime}}}", 25, -10, 10)
  h_outMu_y = ROOT.TH1F("hReal_outMu_y", "; #DeltaY_{#mu^{#prime}} / #sigma_{Y_{#mu^{#prime}}}", 25, -10, 10)
  h_outMu_px = ROOT.TH1F("hReal_outMu_px", "; #Delta(P_{#mu^{#prime}})_{X} / #sigma_(P_{#mu^{#prime}})_{X}", 25, -10, 10)
  h_outMu_py = ROOT.TH1F("hReal_outMu_py", "; #Delta(P_{#mu^{#prime}})_{Y} / #sigma_(P_{#mu^{#prime}})_{Y}", 25, -10, 10)
  h_outMu_pz = ROOT.TH1F("hReal_outMu_pz", "; #Delta(P_{#mu^{#prime}})_{Z} / #sigma_(P_{#mu^{#prime}})_{Z}", 25, -10, 10)

  # Loop over events and fill histograms
  for event in chain: 
    nu = event.nu # using the unfitted version here, honestly idk what to use here 
    if 10 < nu < 32: 
      h_outMu_x.Fill((event.outMu_vec.X() - event.outMuFit_vec.X()) / event.outMu_sigmaX)
      h_outMu_y.Fill((event.outMu_vec.Y() - event.outMuFit_vec.Y()) / event.outMu_sigmaY)
      h_outMu_px.Fill((event.outMu_TL.Px() - event.outMuFit_TL.Px()) / event.outMu_sigmaPx)
      h_outMu_py.Fill((event.outMu_TL.Py() - event.outMuFit_TL.Py()) / event.outMu_sigmaPy)
      h_outMu_pz.Fill((event.outMu_TL.Pz() - event.outMuFit_TL.Pz()) / event.outMu_sigmaPz)

  # Save the histograms to a multi-page PDF, four plots per page
  histograms = (
    h_outMu_x,
    h_outMu_y,
    h_outMu_px,
    h_outMu_py,
    h_outMu_pz,
  )
  saveHistogramsToPdf(histograms, "pull_distributions.pdf")


# **********************************
# Main function 
def main():
  makePullDist()

if __name__ == "__main__":
  main()
