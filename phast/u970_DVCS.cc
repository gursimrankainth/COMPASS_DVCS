#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h" 
#include "PaSetup.h"  
#include "PaAlgo.h"
#include "PaEvent.h" 
#include "PaHodoHelper.h"
#include "G3part.h" 
#include <utility>

#include "ecal_time_cuts.h"
#include "Tools_Camera.hh"
#include "PAMBH_interface_HepGEN.h"

#include "Toololols_KinFitterInterface.hh"
#include "UConn_Tools.h" 

#include "/afs/cern.ch/user/g/gkainth/phastPackages/xcheck_newTiS/tis_range.cc"
#include "/afs/cern.ch/user/g/gkainth/phastPackages/xcheck_newTiS/tis_range.h"

// ************************************************************************** //
// UserEvent for preselecting exclusive photon events	(DVCS)			            //
// In the selection all possible combinations of:					                    //
// Vertices (incoming and outgoing muons), exclusive photon and recoil proton //
// ************************************************************************** //

extern "C" float prob_(float&, int&); 

// Event statistic counter flags  
EventFlags eventFlags(EventFlags::DVCS); // Create an instance of flags and counters 

// Global user selection flags 
bool verbose_mode = false; // Create an instance of verbose_mode
bool leptoMC      = false; // Set to true for LEPTO to remove events with exclusive photoproduction event topology 

//*****************************************************************************
void UserEvent970(PaEvent & e) { // begin event loop

  // Define constants
  static PaHodoHelper* HodoHelper = nullptr; 
  static TiSRange*     tis_range  = nullptr; 
  static PaCamera*     cam_inst   = nullptr;

  //const double M_pi    = G3partMass[8];  // Pion mass 
  //const double M_gamma = G3partMass[1];  // Photon mass
  const double M_mu = G3partMass[5];  // Muon mass
  const double M_p  = G3partMass[14]; // Proton mass 

  // Declare all objects globally
  static TTree* tree(NULL);     // tree for sotring reconstructed data
  static TTree* tree_gen(NULL); // tree for sotring generated data 

  //
  // Variables to be stored into analysis Ntuple
  //
  // (if you want to extend the list:
  // - add variable here,
  // - add in Ntuple definition below 
  // and do not forget to fill the variable in the code
  //
  //*******************************************
  // Shared variables for real and MC data    
  unsigned long long int Evt; // event number - unique evt number (based on run, spill and evt num in spill) 
  static int LastRun = -1; // store the previous run number (used to reintialize Hodohelper, tis_range if there are multiple runs)
  int    Run;          // run number
  int    EvtInSpill;   // event number in spill 
  int    Spill;        // spill number
  double TimeInSpill;  // time in spill  
  double mean_time;    // track meantime 
  float  Chi2;         // Chi2 of the reconstructed vertex
  int    Nprim;        // Number of tracks in primary vertex (-1 in fot found)
  int    Q_beam;       // Charge of the beam 
  int    trig_mask; 

  static BeamFluxParams beamParams; // Create an instance of BeamFluxParams
  static OutMuParams outMuParams;   // Create an instance of OutMuParams 

  static TVector3 SM2center(0, 0, 1825); // Beam charge based on magnetic field 
  static TVector3 SM2field = PaSetup::Ref().MagField(SM2center);

  TVector3 pVtx_vec;     // position vector for the primary vertex (X, Y, Z)
  TVector3 outMu_vec;    // position vector for the scattered muon (use in pull calculations)
  TVector3 posRingA_vec; // hit position of the proton in CAMERA ring A
  TVector3 posRingB_vec; // hit position of the proton in CAMERA ring B

  TLorentzVector inMu_TL;                     // energy-momentum four vector of the beam muon
  TLorentzVector outMu_TL;                    // energy-momentum four vector of the scattered muon
  TLorentzVector gamma_TL;                    // energy-momentum four vector of the real photon 
  TLorentzVector cluster_TL;                  // vector used to store the position and energy information of the cluster
  TLorentzVector p_camera_TL;                 // energy-momentum four vector of the proton measured by camera
  TLorentzVector gammaLow_TL;                 // energy-momentum four vector of the low energy photon (used to find pi0)
  TLorentzVector clusterLow_TL;               // vector used to store the position and energy information of the low-energy cluster
  const TLorentzVector targ_TL(0, 0, 0, M_p); // energy momentum four-vector of the target proton 

  static const int ecal0id = PaSetup::Ref().iCalorim("EC00P1__"); //ECal 0 ID
  static const int ecal1id = PaSetup::Ref().iCalorim("EC01P1__"); // ECal 1 ID
  static const int ecal2id = PaSetup::Ref().iCalorim("EC02P1__"); // ECal 2 ID 
  const float EC0_thr = 4; // ECAL thresholds for exclusive event selection
  const float EC1_thr = 5;
  const float EC2_thr = 10;

  int low_calo  = -1; // calorimeter ID for low energy cluster used for pi0 reconstruction 
  int excl_calo = -1; // calorimeter ID for exclusive photon cluster  

  double y;      // fractional energy loss of the incoming muon 
  double t;      // four-momentum transfer squared of the nucleon 
  double nu;     // energy of the virtual photon 
  double Q2;     // four-momentum transfer sqaured (Q is the four-momentum transferred between the incoming and outgoing lepton)
  double W2;     // effective mass of the final state hadron system squared 
  double xbj;    // measure of the elasticity of the scattering process 
  double phi_gg; // azimuthal angle between the virtual and real photon production planes 

  // Exclusivity variables 
  double delta_phi; 
  double delta_pt;
  double delta_Z;
  double M2x;

  // Kinematic fit;
  TVector3 pVtxFit_vec; 
  TVector3 outMuFit_vec;
  TVector3 posRingAFit_vec; 
  TVector3 posRingBFit_vec; 
  TVector3 clusterFit_vec; // x, y and E of the cluster post kinematic fit 

  TLorentzVector inMuFit_TL; 
  TLorentzVector outMuFit_TL; 
  TLorentzVector targetFit_TL;  
  TLorentzVector gammaFit_TL; 
  TLorentzVector protonFit_TL;  

  double inMu_sigmaX; 
  double inMu_sigmaY;
  double inMu_sigmaPx;
  double inMu_sigmaPy;
  double inMu_sigmaPz;

  double outMu_sigmaX;
  double outMu_sigmaY;
  double outMu_sigmaPx;
  double outMu_sigmaPy;
  double outMu_sigmaPz;

  double gamma_sigmaX;
  double gamma_sigmaY;
  double gamma_sigmaE;

  double proton_sigmaP; 
  double proton_sigmaTheta;
  double proton_sigmaPhi;

  double ringA_sigmaR;
  double ringA_sigmaPhi;
  double ringA_sigmaZ;

  double ringB_sigmaR;
  double ringB_sigmaPhi;
  double ringB_sigmaZ;

  double chi2_fit; // chi2 of the fit  
  int    ndf_fit;  // ndf of the fit 

  double y_fit;
  double nu_fit;
  double Q2_fit; 
  double t_fit;  
  double xbj_fit;
  double phi_gg_fit;  

  //*******************************************
  // HEPGEN BH data (prior to event selection) -> used for acceptance study   
  double phi_gg_gen;
  double weight_all; 
  double weight_DVCS; 
  double weight_BH; 
  double weight_Interference; 
  double phase_fac; 
  double weight_PAMBH; 

  TLorentzVector inMuGen_TL; 
  TLorentzVector outMuGen_TL; 
  TLorentzVector gammaGen_TL; 
  TLorentzVector protonGen_TL;  

  double y_gen;  
  double nu_gen; 
  double Q2_gen;  
  double xbj_gen; 
  double t_gen;  
  double W2_gen; 

  // Generated event selection flag (used later to filter the sample)
  bool genInTarget = false; // generated event vertex is in the target  

  //*******************************************
  // Event selection flags (not the same as statistic counter flags)
  bool trig_MT   = false; // Middle trigger 
  bool trig_LT   = false; // Ladder trigger 
  bool trig_OT   = false; // Outer trigger 
  bool trig_LAST = false;
  bool trig_flag = false; 
  bool TiS_flag  = true;
  bool hodoPass  = false; 
  bool fit_conv  = false; 
  bool save_evt  = false; 
  bool exclLEPTO = false; 

  //*****************************************************************************
  static bool first(true);
  if (first) { // histograms and Ntuples booking block
    
      //
      // Ntuple definition 
      //
      //*******************************************

      tree = new TTree("USR970","User 970 DVCS NTuple"); // name (has to be unique) and title of the Ntuple
      // Basic event information
      tree->Branch("Run", &Run, "Run/I");
      tree->Branch("Evt", &Evt, "Evt/l");
      tree->Branch("Chi2", &Chi2, "Chi2/F");
      tree->Branch("Spill", &Spill, "Spill/I");
      tree->Branch("Q_beam", &Q_beam, "Q_beam/I");
      tree->Branch("TiS_flag", &TiS_flag, "TiS_flag/O");
      tree->Branch("mean_time", &mean_time, "mean_time/D");
      tree->Branch("exclLEPTO", &exclLEPTO, "exclLEPTO/O");
      // Trigger information
      tree->Branch("trig_MT", &trig_MT, "trig_MT/O");
      tree->Branch("trig_LT", &trig_LT, "trig_LT/O");
      tree->Branch("trig_OT", &trig_OT, "trig_OT/O");
      tree->Branch("trig_LAST", &trig_LAST, "trig_LAST/O");
      tree->Branch("trig_flag", &trig_flag, "trig_flag/O");
      tree->Branch("hodoPass", &hodoPass, "hodoPass/O");
      // Particle/vertex vectors 
      tree->Branch("inMu_TL", &inMu_TL);
      tree->Branch("outMu_TL", &outMu_TL);
      tree->Branch("gamma_TL", &gamma_TL);
      tree->Branch("cluster_TL", &cluster_TL);
      tree->Branch("p_camera_TL", &p_camera_TL);
      tree->Branch("gammaLow_TL", &gammaLow_TL);
      tree->Branch("clusterLow_TL", &clusterLow_TL);
      tree->Branch("pVtx_vec", &pVtx_vec); 
      tree->Branch("outMu_vec", &outMu_vec);
      tree->Branch("posRingA_vec", &posRingA_vec);
      tree->Branch("posRingB_vec", &posRingB_vec);
      // Kinematic variables 
      tree->Branch("y", &y, "y/D");
      tree->Branch("t", &t, "t/D");
      tree->Branch("nu", &nu, "nu/D");
      tree->Branch("Q2", &Q2, "Q2/D");
      tree->Branch("W2", &W2, "W2/D");
      tree->Branch("xbj", &xbj, "xbj/D");
      tree->Branch("phi_gg", &phi_gg, "phi_gg/D");
      // Invariant mass of visible pi0
      tree->Branch("low_calo", &low_calo, "low_calo/I");
      tree->Branch("excl_calo", &excl_calo, "excl_calo/I");
      // Exclusivity variables      
      tree->Branch("M2x", &M2x, "M2x/D");
      tree->Branch("delta_Z", &delta_Z, "delta_Z/D");
      tree->Branch("delta_pt", &delta_pt, "delta_pt/D");
      tree->Branch("delta_phi", &delta_phi, "delta_phi/D");
      // Kinematic fit vectors 
      tree->Branch("inMuFit_TL", &inMuFit_TL);
      tree->Branch("outMuFit_TL", &outMuFit_TL);
      tree->Branch("gammaFit_TL", &gammaFit_TL);
      tree->Branch("protonFit_TL", &protonFit_TL);
      tree->Branch("targetFit_TL", &targetFit_TL);
      tree->Branch("pVtxFit_vec", &pVtxFit_vec);
      tree->Branch("clusterFit_vec", &clusterFit_vec);
      tree->Branch("outMuFit_vec", &outMuFit_vec);
      tree->Branch("posRingAFit_vec", &posRingAFit_vec);
      tree->Branch("posRingBFit_vec", &posRingBFit_vec);
      // Kinematic fit 
      tree->Branch("y_fit", &y_fit, "y_fit/D");
      tree->Branch("t_fit", &t_fit, "t_fit/D");
      tree->Branch("nu_fit", &nu_fit, "nu_fit/D");
      tree->Branch("Q2_fit", &Q2_fit, "Q2_fit/D");
      tree->Branch("xbj_fit", &xbj_fit, "xbj_fit/D");
      tree->Branch("ndf_fit", &ndf_fit, "ndf_fit/I");
      tree->Branch("chi2_fit", &chi2_fit, "chi2_fit/D");
      tree->Branch("fit_conv", &fit_conv, "fit_conv/O");
      tree->Branch("phi_gg_fit", &phi_gg_fit, "phi_gg_fit/D");
      tree->Branch("inMu_sigmaX", &inMu_sigmaX, "inMu_sigmaX/D");
      tree->Branch("inMu_sigmaY", &inMu_sigmaY, "inMu_sigmaY/D");
      tree->Branch("inMu_sigmaPx", &inMu_sigmaPx, "inMu_sigmaPx/D");
      tree->Branch("inMu_sigmaPy", &inMu_sigmaPy, "inMu_sigmaPy/D");
      tree->Branch("inMu_sigmaPz", &inMu_sigmaPz, "inMu_sigmaPz/D");
      tree->Branch("outMu_sigmaX", &outMu_sigmaX, "outMu_sigmaX/D");
      tree->Branch("outMu_sigmaY", &outMu_sigmaY, "outMu_sigmaY/D");
      tree->Branch("outMu_sigmaPx", &outMu_sigmaPx, "outMu_sigmaPx/D");
      tree->Branch("outMu_sigmaPy", &outMu_sigmaPy, "outMu_sigmaPy/D");
      tree->Branch("outMu_sigmaPz", &outMu_sigmaPz, "outMu_sigmaPz/D");
      tree->Branch("gamma_sigmaX", &gamma_sigmaX, "gamma_sigmaX/D");
      tree->Branch("gamma_sigmaY", &gamma_sigmaY, "gamma_sigmaY/D");
      tree->Branch("gamma_sigmaE", &gamma_sigmaE, "gamma_sigmaE/D");
      tree->Branch("proton_sigmaP", &proton_sigmaP, "proton_sigmaP/D");
      tree->Branch("proton_sigmaPhi", &proton_sigmaPhi, "proton_sigmaPhi/D");
      tree->Branch("proton_sigmaTheta", &proton_sigmaTheta, "proton_sigmaTheta/D");
      tree->Branch("ringA_sigmaR", &ringA_sigmaR, "ringA_sigmaR/D");
      tree->Branch("ringA_sigmaZ", &ringA_sigmaZ, "ringA_sigmaZ/D");
      tree->Branch("ringA_sigmaPhi", &ringA_sigmaPhi, "ringA_sigmaPhi/D");
      tree->Branch("ringB_sigmaR", &ringB_sigmaR, "ringB_sigmaR/D");
      tree->Branch("ringB_sigmaZ", &ringB_sigmaZ, "ringB_sigmaZ/D");
      tree->Branch("ringB_sigmaPhi", &ringB_sigmaPhi, "ringB_sigmaPhi/D");
      // HEPGEN BH Weights 
      tree->Branch("phase_fac", &phase_fac, "phase_fac/D");
      tree->Branch("weight_BH", &weight_BH, "weight_BH/D");
      tree->Branch("weight_all", &weight_all, "weight_all/D");
      tree->Branch("weight_DVCS", &weight_DVCS, "weight_DVCS/D");
      tree->Branch("weight_PAMBH", &weight_PAMBH, "weight_PAMBH/D");
      tree->Branch("weight_Interference", &weight_Interference, "weight_Interference/D");

      //*******************************************
      tree_gen = new TTree("USR970_gen","User 970 DVCS NTuple (Generated)");
      // Basic event information
      tree_gen->Branch("Run", &Run, "Run/I");
      tree_gen->Branch("Evt", &Evt, "Evt/l");
      tree_gen->Branch("Spill", &Spill, "Spill/I");
      tree_gen->Branch("Q_beam", &Q_beam, "Q_beam/I");
      // Particle/vertex vectors 
      tree_gen->Branch("inMuGen_TL", &inMuGen_TL);
      tree_gen->Branch("outMuGen_TL", &outMuGen_TL);
      tree_gen->Branch("gammaGen_TL", &gammaGen_TL);
      tree_gen->Branch("protonGen_TL", &protonGen_TL);
      // Kinematic variables 
      tree_gen->Branch("y_gen", &y_gen, "y_gen/D");
      tree_gen->Branch("t_gen", &t_gen, "t_gen/D");
      tree_gen->Branch("nu_gen", &nu_gen, "nu_gen/D");
      tree_gen->Branch("Q2_gen", &Q2_gen, "Q2_gen/D");
      tree_gen->Branch("W2_gen", &W2_gen, "W2_gen/D");
      tree_gen->Branch("xbj_gen", &xbj_gen, "xbj_gen/D");
      tree_gen->Branch("phi_gg_gen", &phi_gg_gen, "phi_gg_gen/D");
      // HEPGEN BH Weights 
      tree_gen->Branch("phase_fac", &phase_fac, "phase_fac/D");
      tree_gen->Branch("weight_BH", &weight_BH, "weight_BH/D");
      tree_gen->Branch("weight_all", &weight_all, "weight_all/D");
      tree_gen->Branch("weight_DVCS", &weight_DVCS, "weight_DVCS/D");
      tree_gen->Branch("weight_PAMBH", &weight_PAMBH, "weight_PAMBH/D");
      tree_gen->Branch("weight_Interference", &weight_Interference, "weight_Interference/D");
      // Flags 
      tree_gen->Branch("genInTarget", &genInTarget, "genInTarget/O");

      first = false;
  } // end of histogram booking

  //*****************************************************************************
  // Assign names to trigger bits
  enum trigger {   
    Tiger = 1<<0,
    MT = 1<<1,
    LT = 1<<2,
    OT = 1<<3,
    CT = 1<<4,
    IV = 1<<5,
    HaloT = 1<<6,
    BT = 1<<7,
    Tiger_only = 1<<8,
    LAST = 1<<9,
    TRand = 1<<10,
    NRand = 1<<11
  };

  trig_mask = e.TrigMask(); 
  trig_mask = trig_mask&2047;

  std::string trigCheck = ""; // empty string to store trigger information for debugging 
  if (trig_mask & MT) {
      trigCheck += "MT ";
      trig_flag = true;
      trig_MT = true; 
  }
  if (trig_mask & LT) {
      trigCheck += "LT ";
      trig_flag = true;
      trig_LT = true;
  }
  if (trig_mask & OT) {
      trigCheck += "OT ";
      trig_flag = true;
      trig_OT = true;
  }
  if (trig_mask & LAST) {
      trigCheck += "LAST ";
      trig_flag = true; 
      trig_LAST = true;
  }

  //*******************************************
  // Initialize variables, extra flags and check time in spill (time in spill cut is applied later not here)   
  // Check flags are only used to check statistics not for final event selection 
  eventFlags.createFlag("singleTrack_flag", "No. of events where primary vertex only has one outgoing track");
  eventFlags.createFlag("Q2_DIS_flag", "No. of events where Q2 > 0.5 (DIS)");
  //eventFlags.createFlag("y_DIS_flag", "No. of events where 0.01 < y < 0.99");

  eventFlags.createFlag("clAll_flag", "No. of events that have any clusters");
  eventFlags.createFlag("clNeutral_flag", "No. of events where clusters are not associated with charged tracks");
  eventFlags.createFlag("clTime_flag", "No. of events where cluster timing is within requirements");
  eventFlags.createFlag("lowECl_flag", "No. of events with low energy photons");
  eventFlags.createFlag("nCls_flag", "No. of events where high energy clusters are found in ECal 0, 1 or 2 only");
  eventFlags.createFlag("singleCl_flag", "No. of events where there is only a single high energy cluster in the ECals");

  eventFlags.createFlag("protonAll_flag", "No. of events with proton candidates");
  eventFlags.createFlag("proton_flag", "No. of events where proton candidates have 0.1 < beta < 1");
  eventFlags.createFlag("delta_pt_flag", "No. of events where |delta_pt| < 0.3 GeV/c");
  eventFlags.createFlag("delta_phi_flag", "No. of events where |delta_phi| < 0.4 rad");
  eventFlags.createFlag("delta_Z_flag", "No. of events where |delta_Z| < 16 cm");
  eventFlags.createFlag("M2x_flag", "No. of events where |(M_x)^2| < 0.3 (GeV/c^2)^2");

  eventFlags.createFlag("Q2Fit_flag", "No. of events where 1 < Q2_fit < 10");
  eventFlags.createFlag("yFit_flag", "No. of events where 0.05 < y_fit < 0.95"); //final DVCS cut is at 0.9 not 0.95 
  eventFlags.createFlag("tFit_flag", "No. of events where -0.08 < t_fit < -0.64");
  eventFlags.createFlag("nuFit_flag", "No. of events where 10 < nu_fit < 144");
  eventFlags.createFlag("kinFitAll_flag", "No. of events surviving all kinematic cuts (Q2, y, t, nu)");
  
  eventFlags.createFlag("nExclCombo_flag", "No. of events where at least one vertex, photon and proton combination satisfies 4/5 exclusivity cuts");
  eventFlags.createFlag("nExclComboPi0_flag", "No. of events with pi0 candidates which satisfy at least 4/5 exclsuvity cuts");
  eventFlags.createFlag("nComboPi0_flag", "No. of events (pi0 MC) with pi0 candidates - no pi0 exclusivity cuts"); 
  eventFlags.createFlag("nPi0Save_flag", "No. of events (pi0 MC) with pi0 candidates which satisfy at least 4/5 exclsuvity cuts"); 
  
  eventFlags.createFlag("allGen_flag", "No. of all generated events");
  eventFlags.createFlag("saveGen_flag", "No. of generated events which are saved"); 

  eventFlags.resetFlags(); // Reset all event statistic counter flags to false

  Run         = e.RunNum();
  Evt         = e.UniqueEvNum();
  EvtInSpill  = e.EvInSpill(); 
  Spill       = e.SpillNum(); 
  TimeInSpill = e.TimeInSpill();  

  if (Run != LastRun) { // Reinitialize HodoHelper and tis_range only if the run number changes 
      HodoHelper = & PaHodoHelper::Init("/afs/cern.ch/user/g/gkainth/phast/dat/trigger_config/2016", true);  
      //OLD: tis_range  = new TiSRange("/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files");
      tis_range  = new TiSRange("/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files_slot7.1");
      LastRun = Run;  // Update LastRun to the current run number
  }

  if (!e.IsMC()) {
    TiS_flag = tis_range->CheckWindow(Run, Spill, TimeInSpill); // check the time in spill for real data
  }

  cam_inst = & PaCamera::GetInstance();
  //cam_inst->NewEvent(e, e.IsMC()); // only have it like this for the x-check for the full analysis use cam_inst->NewEvent(e); 
  cam_inst->NewEvent(e); 

  //*******************************************
  // Reset all generated variables FIRST
  resetGenVars(
    y_gen, Q2_gen, nu_gen, xbj_gen, t_gen, W2_gen, phi_gg_gen,
    weight_all, weight_DVCS, weight_BH, weight_Interference,
    weight_PAMBH, phase_fac,
    genInTarget
  );

  //*******************************************
  // Get MC weights (for HEPGEN BH)
  bool MCevent = false;
  if (e.IsMC()) { // Begin loop over MC data 
    MCevent = true;
    eventFlags.setFlagByName("allGen_flag", true);
    NLUDATA ld; 
    if (e.MCgen(ld)) { // Begin loop to get event weights 
      weight_all  = ld.uservar[2]; 
      weight_DVCS = ld.uservar[15];  
      weight_BH   = ld.uservar[16];   
      phase_fac   = ld.uservar[9]; 
    } // End loop to get event weights 
  } // End loop over MC data 

  // Check for exclusive event topology in the case of LEPTO
  enum class LujetCheckStatus {NOT_MC, I};
  if (leptoMC && exclLepto(e)) {
    exclLEPTO = true;  // TRUE = exclusive pi0
  }
  
  //*******************************************  
  // Loop over reconstructed vertices in the event - REAL and MC RECONSTRUCTED DATA 
  eventFlags.setFlagByName("allEvts_flag", true); 
  bool genFilled = false;
  for (int iv = 0; iv < e.NVertex(); iv++) { // Begin loop over vertices
    //******************************************* 
    // Store info about primary vertex (if found) 
    const PaVertex &v = e.vVertex(iv); 
    if (!v.IsPrimary()) continue; // skip any vertices that are not primary 
    eventFlags.setFlagByName("pVtx_flag", true);
    pVtx_vec.SetXYZ(v.Pos(0), v.Pos(1), v.Pos(2));
    Chi2  = v.Chi2(); 
    Nprim = v.NOutParticles(); // number of tracks in vertex

    //*******************************************
    // Store info about incoming muon beam (inMu)
    Q_beam = SM2field(1) < 0 ? 1 : -1;

    // Check that there is an incoming muon associated with the vertex
    int i_beam = v.InParticle();  
    if (i_beam == -1) continue;

    // Check that the beam has a track associated with it
    PaParticle beam = e.vParticle(i_beam);
    int it_beam = beam.iTrack();
    if (it_beam == -1) continue;

    // Check that the track has parameters
    PaTrack beam_track = e.vTrack(it_beam);
    if (beam_track.NTPar() == 0) continue;
    eventFlags.setFlagByName("inMuTrack_flag", true);

    // Check that the beam was first measured before the target
    if (beam_track.ZFirst() >= -78.5) continue;
    eventFlags.setFlagByName("zFirst_flag", true);

    // Check that the beam momentum falls within acceptable range
    double inMu_mom = beam_track.vTPar(0).Mom();
    if (inMu_mom < 140.0 || inMu_mom > 180.0) continue;
    eventFlags.setFlagByName("momRange_flag", true);

    // Check that the beam momentum error falls within acceptable range
    double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
    if (inMu_momErr > 0.025*inMu_mom) continue;
    eventFlags.setFlagByName("momErr_flag", true); 

    // Check that the beam is detected by detectors along the beamline
    int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
    int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
    int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 
    
    if (!e.IsMC()) { // skip BMS check for MC events 
      if (nhits_BMS < 3) continue;
      eventFlags.setFlagByName("BMS_flag", true);
    } else {eventFlags.setFlagByName("BMS_flag", true);}

    if (nhits_FI < 2) continue;
    eventFlags.setFlagByName("FI_flag", true);

    if (nhits_SI < 3) continue;
    eventFlags.setFlagByName("SI_flag", true);

    // Check that the beam crosses the full target length 
    PaTPar par_beam = beam.ParInVtx(iv); // beam parameters at the vertex 
    if (!PaAlgo::CrossCells(beam_track.vTPar(0), Run, beamParams.Rmax, beamParams.Ymax, beamParams.tgt_zmin, beamParams.tgt_zmax, beamParams.RmaxMC)) continue;
    eventFlags.setFlagByName("crossCells_flag", true);

    // Check that the track meantime is within flux requirements
    mean_time = beam_track.MeanTime(); 
    if (!e.IsMC()) { // ignore meantime check for MC????? 
      if (std::fabs(mean_time) >= 2) continue;
      eventFlags.setFlagByName("meantime_flag", true);
    } else {eventFlags.setFlagByName("meantime_flag", true);}

    //*******************************************
    // Save to generated tree 
    if (!genFilled && e.IsMC() && !leptoMC) { // Begin loop over MC data 
      for (int imcv = 0; imcv < e.NMCvertex(); imcv++) { 
        const PaMCvertex &vmc = e.vMCvertex(imcv);
        if (!vmc.IsPrimary()) continue; 

        if (vmc.NMCtrack() != 4) continue; 

        const PaMCtrack &t_beam_mc   = e.vMCtrack(vmc.iMCtrack(0));
        const PaMCtrack &t_outMu_mc  = e.vMCtrack(vmc.iMCtrack(1));
        const PaMCtrack &t_gamma_mc  = e.vMCtrack(vmc.iMCtrack(2));
        const PaMCtrack &t_proton_mc = e.vMCtrack(vmc.iMCtrack(3));

        const PaTPar &par_beam_mc   = t_beam_mc.ParInVtx();
        const PaTPar &par_outMu_mc  = t_outMu_mc.ParInVtx();
        const PaTPar &par_gamma_mc  = t_gamma_mc.ParInVtx();
        const PaTPar &par_proton_mc = t_proton_mc.ParInVtx();

        inMuGen_TL   = par_beam_mc.LzVec(M_mu);
        outMuGen_TL  = par_outMu_mc.LzVec(M_mu);
        gammaGen_TL  = par_gamma_mc.LzVec(0);
        protonGen_TL = par_proton_mc.LzVec(M_p);

        // Set flag if the gnerated vertex is in the target  
        genInTarget = PaAlgo::InTarget(par_beam,'O', Run, beamParams.Rmax, beamParams.Ymax, beamParams.tgt_zmin, beamParams.tgt_zmax, beamParams.RmaxMC); 

        t_gen   = (targ_TL - protonGen_TL) * (targ_TL - protonGen_TL);
        y_gen   = (inMuGen_TL.E() - outMuGen_TL.E()) / inMuGen_TL.E(); 
        Q2_gen  = PaAlgo::Q2 (inMuGen_TL, outMuGen_TL);
        nu_gen  = inMuGen_TL.E() - outMuGen_TL.E();
        W2_gen  = PaAlgo::W2 (inMuGen_TL, outMuGen_TL);
        xbj_gen = PaAlgo::xbj (inMuGen_TL, outMuGen_TL);

        double E_beam_gen = inMuGen_TL.E(); 
        phi_gg_gen = phiRV(inMuGen_TL, outMuGen_TL, protonGen_TL, gammaGen_TL, true);
        weight_PAMBH = Weight_Pam_BH(xbj_gen, Q2_gen, phi_gg_gen, t_gen, E_beam_gen, phase_fac); 

        tree_gen->Fill();
        eventFlags.setFlagByName("saveGen_flag", true);
        genFilled = true;
        break; // ensure that event only fills once per vertex
      }
    } // End loop over MC data

    // TiS flag is used only to check DVCS statistics -> cut is in postPhast script 
    if (TiS_flag) {  
      eventFlags.setFlagByName("timeInSpill_flag", true); 
    }  
    
    //*******************************************
    // Store info about scattered muon (outMu)
    static PaParticle outMu; 
    static PaTrack outMu_track; 
    static PaTPar par_outMu; 

    // outMu_flag will be true as long as any muon is found ( and all other conditions satisfied), even if it doesnt pass the full Hodoscope check 
    bool outMu_flag = outMuCheck(e, v, iv, Run, beam, HodoHelper, trig_flag, TiS_flag, outMuParams, 
                                outMu, outMu_track, par_outMu, eventFlags); 
    if (!outMu_flag) continue; // ignore events that dont satisfy requirements for scattered muon

    outMu_vec.SetXYZ(par_outMu.X(), par_outMu.Y(), par_outMu.Z());

    int i_omu_check_hodo = HodoHelper->iMuPrim(v, false, false, true, true, 15, true, true);
    if (i_omu_check_hodo != -1) {hodoPass = true;}

    //*******************************************
    // Kinematic variables ... (1/2)
    inMu_TL  = par_beam.LzVec(M_mu); 
    outMu_TL = par_outMu.LzVec(M_mu);  

    Q2  = PaAlgo::Q2 (inMu_TL, outMu_TL); //Q2 = -(inMu_TL - outMu_TL).M2();
    y   = (inMu_TL.E() - outMu_TL.E()) / inMu_TL.E();
    nu  = (inMu_TL.E() - outMu_TL.E());
    W2  = PaAlgo::W2 (inMu_TL, outMu_TL);
    xbj = PaAlgo::xbj (inMu_TL, outMu_TL); //xbj = Q2/(2*q*targ_TL);

    //*******************************************
    // Exclusive selection starts here  ...  
    // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
    if (Nprim != 1) continue; 
    eventFlags.setFlagByName("singleTrack_flag", true);

    // Current kinematic cuts will be tightened after the kinematically constrained fit is applied 
    if (Q2 < 0.5) continue; // inclusive Q2 cut
    eventFlags.setFlagByName("Q2_DIS_flag", true);   

    //if (y < 0.01 || y > 0.99) continue; // inclusive y cut
    //eventFlags.setFlagByName("y_DIS_flag", true);

    //*******************************************
    // Store info about real photon - check ECals for single neutral cluster (can be from any exclusive event (DVCS or BH))
    // ! CURRENTLY LOOKING FOR ANY NEUTRAL CLUSTERS, NOT SPECIFICALLY FOR ONE !
    std::vector<int> cl_id; // Store cluster index for valid exclusive clusters
    std::vector<int> pi0_cl_id; // Store cluster index for low-energy clusters 

    for (int iclus = 0; iclus < e.NCaloClus(); iclus++) { // Begin loop over clusters 
      eventFlags.setFlagByName("clAll_flag", true);
      const PaCaloClus & cl = e.vCaloClus(iclus);
      int icalo = cl.iCalorim();

      if (cl.iTrack() != -1) continue; // Require neutral
      eventFlags.setFlagByName("clNeutral_flag", true);

      if (!EcalTimeCut(beam_track, cl, Run) && !e.IsMC()) continue;
      eventFlags.setFlagByName("clTime_flag", true);

      // Remove any clusters which do not have the required ECal ID 
      if ((icalo != ecal0id) && (icalo != ecal1id) && (icalo != ecal2id)) continue;

      float E = cl.E();
      if ((icalo == ecal0id && E >= EC0_thr) ||
          (icalo == ecal1id && E >= EC1_thr) ||
          (icalo == ecal2id && E >= EC2_thr)) {
        cl_id.push_back(iclus);
      } else {
        pi0_cl_id.push_back(iclus);
      }

    } // End loop over clusters 

    if (cl_id.size() != 0) {
      eventFlags.setFlagByName("nCls_flag", true);
    }

    if (pi0_cl_id.size() != 0) {
      eventFlags.setFlagByName("lowECl_flag", true);
    }

    // ! NOW CHECK THAT THERE IS ONLY ONE NEUTRAL CLUSTER !
    if (cl_id.size() != 1) continue; // skip event if there is more than one exclusive cluster
    eventFlags.setFlagByName("singleCl_flag", true);

    // Calculate TLorentz vectors for single high energy photon
    buildClusterVecs(e, v, cl_id[0], gamma_TL, cluster_TL);
    excl_calo = e.vCaloClus(cl_id[0]).iCalorim(); // get the ECal ID for the exclusive cluster

    //*******************************************
    // Store info about scattered proton candidates (check CAMERA)
    // ! CURRENTLY LOOKING FOR ANY PROTONS, NOT SPECIFICALLY FOR ONE !
    vector <CameraProton> proton_candidates = cam_inst->GetGoodCandidates(v); // vector holding all proton candidates

    TVector3 R_vtx;
    R_vtx.SetXYZ(v.Pos(0),v.Pos(1),v.Pos(2));

    // For exclusive photoproduction selection - real data + HEPGEN BH
    bool betaPassed = false; 
    TLorentzVector p_spec_TL = targ_TL + inMu_TL - outMu_TL - gamma_TL;
    double pt_spec = p_spec_TL.Pt(); 
    double phi_spec = p_spec_TL.Phi();

    if (proton_candidates.size() == 0) continue;
    for (auto proton: proton_candidates) { // Begin loop over proton candidates
      eventFlags.setFlagByName("protonAll_flag", true); 

      bool DeltaPtPassed  = false; 
      bool DeltaPhiPassed = false; 
      bool DeltaZPassed   = false; 
      bool DeltaM2xPassed = false; 
      bool tPassed        = false; 

      double beta = proton.beta;
      if (beta >= 0.1 && beta <= 1) {betaPassed = true;}
      if (!betaPassed) continue;
      eventFlags.setFlagByName("proton_flag", true);

      p_camera_TL = proton.p4; 
      if (p_camera_TL.Mag() == 0) continue; // ignore events where there is no TL vector

      t = (targ_TL - p_camera_TL) * (targ_TL - p_camera_TL); 
      phi_gg = phiRV(inMu_TL, outMu_TL, p_camera_TL, gamma_TL, true);

      //*******************************************
      // Check that all combintations of the vertex, photon and proton satisfy exclusivity conditions 
      // 4 exclusivity variables: detla_phi, delta_pt, delta_Z, M2x
      double pt_camera = p_camera_TL.Pt(); 
      delta_pt = pt_camera - pt_spec; // transverse momentum
      double phi_camera = p_camera_TL.Phi();
      delta_phi = phi_camera - phi_spec; // azimuthal angle

      posRingA_vec = proton.Ahit.vec;
      posRingB_vec = proton.Bhit.vec;
      double Z_inter = cam_inst->GetZA(R_vtx, posRingA_vec, posRingB_vec, phi_camera);
      delta_Z = posRingA_vec.Z() - Z_inter; // z position of the hits in the inner CAMERA ring

      M2x = (p_spec_TL - p_camera_TL) * (p_spec_TL - p_camera_TL);

      int protonPassCount = 0;
      if ((DeltaPtPassed  = std::fabs(delta_pt) <= 0.3)) protonPassCount++;
      if ((DeltaZPassed   = std::fabs(delta_Z) <= 16)) protonPassCount++;
      if ((DeltaM2xPassed = std::fabs(M2x) <= 0.3)) protonPassCount++;
      if ((tPassed = (t < -0.08 && t > -0.64))) protonPassCount++;

      if (delta_phi >= -0.4 && delta_phi <= 0.4) {
          DeltaPhiPassed = true;
          protonPassCount++;
      } else if ((delta_phi + 2 * TMath::Pi()) >= -0.4 && (delta_phi + 2 * TMath::Pi()) <= 0.4) {
          delta_phi += 2 * TMath::Pi();
          DeltaPhiPassed = true;
          protonPassCount++;
      } else if ((delta_phi - 2 * TMath::Pi()) >= -0.4 && (delta_phi - 2 * TMath::Pi()) <= 0.4) {
          delta_phi -= 2 * TMath::Pi();
          DeltaPhiPassed = true;
          protonPassCount++;
      }

      if (TiS_flag && (eventFlags.getFlag("passHodo_flag"))) {
        if (Q2 >= 1 && Q2 <= 10) {
          eventFlags.setFlagByName("Q2_DVCS_flag", true);
          if (y >= 0.05 && y <= 0.9) {
            eventFlags.setFlagByName("y_DVCS_flag", true);
            if (DeltaPtPassed) {
              eventFlags.setFlagByName("delta_pt_flag", true); 
              if (DeltaPhiPassed) {
                eventFlags.setFlagByName("delta_phi_flag", true);
                if (DeltaZPassed) {
                  eventFlags.setFlagByName("delta_Z_flag", true);
                  if (DeltaM2xPassed) {
                    eventFlags.setFlagByName("M2x_flag", true);
                  }
                }
              }
            }
          }
        }
      }

      //*******************************************
      // Perform the kinematic fit for current combination and save results 
      static Fitter* FitInterface = &(Fitter::GetInstance());
      FitInterface->Init(R_vtx, beam_track, outMu_track, p_camera_TL, posRingA_vec, posRingB_vec); 
      FitInterface->Add_Photon(e.vCaloClus(cl_id[0]));
      FitInterface->SetupFit();
      FitInterface->DoFit(0, 1000);

      pVtxFit_vec     = *(FitInterface->GetVertex()->getCurr3Vec());
      posRingAFit_vec = *(FitInterface->GetHitA()->getCurr3Vec());
      posRingBFit_vec = *(FitInterface->GetHitB()->getCurr3Vec()); 
      clusterFit_vec  = *(FitInterface->GetOutPhotons()[0]->getCurr3Vec());
      outMuFit_vec  = *(FitInterface->GetMuonOut()->getCurr3Vec());

      inMuFit_TL   = *(FitInterface->GetMuonIn()->getCurr4Vec());
      outMuFit_TL  = *(FitInterface->GetMuonOut()->getCurr4Vec());
      protonFit_TL = *(FitInterface->GetProtonOut()->getCurr4Vec());
      targetFit_TL = *(FitInterface->GetProtonTarget()->getCurr4Vec());
      gammaFit_TL  = *(FitInterface->GetOutPhotons()[0]->getCurr4Vec());

      inMu_sigmaX  = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(0, 0));
      inMu_sigmaY  = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(1, 1));
      inMu_sigmaPx = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(2, 2));
      inMu_sigmaPy = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(3, 3));
      inMu_sigmaPz = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(4, 4));

      outMu_sigmaX  = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(0, 0));
      outMu_sigmaY  = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(1, 1));
      outMu_sigmaPx = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(2, 2));
      outMu_sigmaPy = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(3, 3));
      outMu_sigmaPz = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(4, 4));

      gamma_sigmaX = TMath::Sqrt((*(FitInterface->GetOutPhotons()[0]->getCovMatrixDeltaY()))(0, 0));
      gamma_sigmaY = TMath::Sqrt((*(FitInterface->GetOutPhotons()[0]->getCovMatrixDeltaY()))(1, 1));
      gamma_sigmaE = TMath::Sqrt((*(FitInterface->GetOutPhotons()[0]->getCovMatrixDeltaY()))(2, 2));

      proton_sigmaP     = TMath::Sqrt((*(FitInterface->GetProtonOut()->getCovMatrixDeltaY()))(0, 0));
      proton_sigmaTheta = TMath::Sqrt((*(FitInterface->GetProtonOut()->getCovMatrixDeltaY()))(1, 1));
      proton_sigmaPhi   = TMath::Sqrt((*(FitInterface->GetProtonOut()->getCovMatrixDeltaY()))(2, 2));

      ringA_sigmaR   = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(0, 0));
      ringA_sigmaPhi = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(1, 1));
      ringA_sigmaZ   = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(2, 2));

      ringB_sigmaR   = TMath::Sqrt((*(FitInterface->GetHitB()->getCovMatrixDeltaY()))(0, 0));
      ringB_sigmaPhi = TMath::Sqrt((*(FitInterface->GetHitB()->getCovMatrixDeltaY()))(1, 1));
      ringB_sigmaZ   = TMath::Sqrt((*(FitInterface->GetHitB()->getCovMatrixDeltaY()))(2, 2));
      
      fit_conv = FitInterface->GetFitOutput(chi2_fit, ndf_fit);  

      Q2_fit = PaAlgo::Q2 (inMuFit_TL, outMuFit_TL); 
      xbj_fit = PaAlgo::xbj (inMuFit_TL, outMuFit_TL);
      y_fit  = (inMuFit_TL.E() - outMuFit_TL.E()) / inMuFit_TL.E(); 
      nu_fit = (inMuFit_TL.E() - outMuFit_TL.E()); 
      t_fit  = (targetFit_TL - protonFit_TL) * (targetFit_TL - protonFit_TL);  
      phi_gg_fit = phiRV(inMuFit_TL, outMuFit_TL, protonFit_TL, gammaFit_TL, MCevent);

      bool Q2_cut = false; 
      bool y_cut  = false; 
      bool t_cut  = false; 
      bool nu_cut = false; 

      if ((Q2_fit > 1 && Q2_fit < 10) || std::isnan(Q2_fit)) {
        Q2_cut = true;
        eventFlags.setFlagByName("Q2Fit_flag", true);
      }
      if ((y_fit > 0.05 && y_fit < 0.95) || std::isnan(y_fit)) {
        y_cut = true; 
        eventFlags.setFlagByName("yFit_flag", true);
      }
      if ((t_fit < -0.08 && t_fit > -0.64) || std::isnan(t_fit)) {
        t_cut = true;
        eventFlags.setFlagByName("tFit_flag", true);
      }
      if ((nu_fit > 10 && nu_fit < 144) || std::isnan(nu_fit)) { 
        nu_cut = true; 
        eventFlags.setFlagByName("nuFit_flag", true);
      }
      if (Q2_cut && y_cut && t_cut && nu_cut) {
        eventFlags.setFlagByName("kinFitAll_flag", true);
      }

      //*******************************************
      // Choose exclusive candidates to save
      if (protonPassCount >= 4) { // Begin loop over candidates which pass 4/5 exclusivity cuts 
        eventFlags.setFlagByName("nExclCombo_flag", true);
        save_evt = true;

        if (!pi0_cl_id.empty()) { // Begin loop over low-energy clusters
          for (auto iLow = std::size_t{0}; iLow < pi0_cl_id.size(); ++iLow) {  
            eventFlags.setFlagByName("nExclComboPi0_flag", true);
            const auto& cl_LowE = e.vCaloClus(pi0_cl_id[iLow]);
            low_calo = cl_LowE.iCalorim();
            buildClusterVecs(e, v, pi0_cl_id[iLow], gammaLow_TL, clusterLow_TL);
            tree->Fill();
          }
        } else {
          tree->Fill();
        } // End loop over low-energy clusters

      } // End loop over candidates which pass 4/5 exclusivity cuts 
    } // End loop over proton candidates

    //*******************************************
    // Debug statements ...
    printDebug("     ");
    printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
    printDebug("    Vertex: (" + std::to_string(pVtx_vec.X()) + ", " + std::to_string(pVtx_vec.Y()) + ", " + std::to_string(pVtx_vec.Z()) + ")");
    printDebug("    mu: P: " + std::to_string(inMu_TL.P()) + " GeV/c, Charge: " + std::to_string(beam.Q()));
    printDebug("    mu': P: " + std::to_string(outMu_TL.P()) + " GeV/c, Charge: " + std::to_string(outMu.Q()));
    printDebug("    Kinematics: Q2: " + std::to_string(Q2) +  " GeV2, y: " + std::to_string(y) + ", W2: " + std::to_string(W2) + " GeV2, x: " + std::to_string(xbj));

  } // End loop over vertices 

  // Increment all counters whose flags are "true"
  eventFlags.incrementCounters(); 
  //Save the event
  if (save_evt) {e.TagToSave();}    

} // End event loop 

void UserJobEnd970() {
  eventFlags.printFlags(); // Print to output stream 
}
