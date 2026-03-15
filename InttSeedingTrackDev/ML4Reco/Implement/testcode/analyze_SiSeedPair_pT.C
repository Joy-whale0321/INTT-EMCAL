#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TNtuple.h>

#include <TLorentzVector.h>

#include <cmath>
#include <string>
#include <vector>

#include "PtCalculator.h"  // SiCaloPt::PtCalculator & friends
R__LOAD_LIBRARY(libPtCalc.so)


const float electron_mass = 0.000511; // GeV/c^2
const float pion_mass     = 0.1396; // GeV/c^2

SiCaloPt::PtCalculator pTCalc;

// ---- Weights(onnx) and Scalers(json) Path ---------------------------
struct DemoPaths
{
  std::string file_dir             = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/"; 
  std::string emd_onnx             = file_dir + "model_MLEMD.onnx"; 
  std::string emd_scaler_json      = "";
  
  std::string eproj_onnx           = file_dir + "model_MLEproj.onnx"; 
  std::string eproj_scaler_json    = file_dir + "scaler_MLEproj.json"; 
  
  std::string combined_onnx        = file_dir + "model_MLCombined.onnx"; 
  std::string combined_scaler_json = "";

  void print(){
    std::cout<<"emd_onnx             : "<<emd_onnx            <<std::endl;
    std::cout<<"emd_scaler_json      : "<<emd_scaler_json     <<std::endl;
    std::cout<<"eproj_onnx           : "<<eproj_onnx          <<std::endl;
    std::cout<<"eproj_scaler_json    : "<<eproj_scaler_json   <<std::endl;
    std::cout<<"combined_onnx        : "<<combined_onnx       <<std::endl;
    std::cout<<"combined_scaler_json : "<<combined_scaler_json<<std::endl;
  };
};

// ---- turn string into optional<string> for Config ------------------------
template<typename Opt>
Opt make_opt(const std::string& s) 
{
  if (s.empty()) return std::nullopt;
  return s;
}

struct PtInput {
  float SiIn[3];
  float SiOut[3];
  float Calo[3];
  float CaloE;
};

float calcPt(PtInput& in)
{
  // for EMD and ML-EMD
  float phi_si   = atan2(in.SiOut[1]-in.SiIn[1], in.SiOut[0]-in.SiIn[0]);
  float phi_calo = atan2(in.Calo[1]-in.SiOut[1], in.Calo[0]-in.SiOut[0]);
  float dphi = phi_calo - phi_si;

  // for EMD
  SiCaloPt::InputEMD inEMD;
  inEMD.EMD_Angle  = dphi; // delta_phi - EM Deflection angle in rad
  inEMD.EMD_Eta    = 0.00; // EMCal Cluster eta
  inEMD.EMD_Radius = 93.5; // EMCal Cluster radius in cm

  auto r = pTCalc.ComputePt(SiCaloPt::Method::MethodEMD, SiCaloPt::AnyInput{inEMD});
  std::cout << "[EMD-analytic] ok=" << r.ok
            << "  pt=" << r.pt_reco
            << "  err=\"" << r.err << "\"\n";

  float pT = r.pt_reco;

  // for ML-EMD
  // 2-d input features:{ dphi_EMD, eta_track }
  std::vector<float> featsMLEMD = {dphi, 0};
  SiCaloPt::InputMLEMD inMLEMD{featsMLEMD};

  r = pTCalc.ComputePt(SiCaloPt::Method::MethodMLEMD, SiCaloPt::AnyInput{inMLEMD});
  std::cout << "[MLEMD-2D] ok=" << r.ok
          << "  pt=" << r.pt_reco
          << "  err=\"" << r.err << "\"\n";


  // for ML-Eproj
  float r_siin  = sqrt(pow(in.SiIn[0], 2) + pow(in.SiIn[1], 2));
  float r_siout = sqrt(pow(in.SiOut[0], 2)+ pow(in.SiOut[1], 2));
  float r_calo  = sqrt(pow(in.Calo[0], 2)+ pow(in.Calo[1], 2));


  // 7-d input features:{ INTT 3/4 layer R, INTT 3/4 layer Z, INTT 5/4 layer R, INTT 5/4 layer Z, Calo R, Calo Z, Calo Energy }
  std::vector<float> featsMLEproj = { r_siin,  in.SiIn[2],   // INTT 3/4 layer
                                      r_siout, in.SiOut[2],  // INTT 5/4 layer
                                      r_calo,  in.Calo[2], in.CaloE}; // Calo R,Z,Energy
  
  SiCaloPt::InputMLEproj inMLEproj{featsMLEproj};
  r = pTCalc.ComputePt(SiCaloPt::Method::MethodMLEproj, SiCaloPt::AnyInput{inMLEproj});
  std::cout << "[MLEproj-7D] ok=" << r.ok
          << "  pt=" << r.pt_reco
          << "  err=\"" << r.err << "\"\n";



  return pT;
}

float demoPt(){

  // ======================== EMD Formula ==========================
  {
    SiCaloPt::InputEMD in;
    in.EMD_Angle  = 0.025;  // delta_phi - EM Deflection angle in rad
    in.EMD_Eta    = 0.00;   // EMCal Cluster eta
    in.EMD_Radius = 93.5;   // EMCal Cluster radius in cm
    
    pTCalc.setParCeta(0.2);   // set C_eta parameter if needed
    pTCalc.setParPower(-1.0); // set Power parameter if needed
    
    auto r = pTCalc.ComputePt(SiCaloPt::Method::MethodEMD, SiCaloPt::AnyInput{in});
    std::cout << "[EMD-analytic] ok=" << r.ok
              << "  pt=" << r.pt_reco
              << "  err=\"" << r.err << "\"\n";

    in.EMD_Angle  = -0.025;  // delta_phi - EM Deflection angle in rad
    r = pTCalc.ComputePt(SiCaloPt::Method::MethodEMD, SiCaloPt::AnyInput{in});
    std::cout << "[EMD-analytic] ok=" << r.ok
              << "  pt=" << r.pt_reco
              << "  err=\"" << r.err << "\"\n";
  }
  
  // ======================== Eproj Formula ==========================
  {
    SiCaloPt::InputEproj in;
    in.Energy_Calo   = 1.8;   // EMCal Cluster energy in GeV
    in.Radius_Calo   = 93.5;  // EMCal Cluster radius in cm
    in.Z_Calo        = 0.0;   // EMCal Cluster z in cm
    in.Radius_vertex = 0.0;   // Vertex radius in cm
    in.Z_vertex      = 0.0;   // Vertex z in cm
    
    auto r = pTCalc.ComputePt(SiCaloPt::Method::MethodEproj, SiCaloPt::AnyInput{in});
    std::cout << "[Eproj-analytic] ok=" << r.ok
      << "  pt=" << r.pt_reco
      << "  err=\"" << r.err << "\"\n";
  }


  // ============= ML：MLEMD (2-d input: 1/dphi_EMD, eta_track = 0 now) ===============
  {
    // 2-d input features:{ dphi_EMD, eta_track }
    std::vector<float> featsMLEMD = {15, 0};
  
    SiCaloPt::InputMLEMD in{featsMLEMD};
    auto r = pTCalc.ComputePt(SiCaloPt::Method::MethodMLEMD, SiCaloPt::AnyInput{in});
    std::cout << "[MLEMD-2D] ok=" << r.ok
      << "  pt=" << r.pt_reco
      << "  err=\"" << r.err << "\"\n";
  }

  // ====== ML：MLEproj (7-d input: INTT 3/4 layer R,Z, INTT 5/4 layer R,Z, Calo R,Z,Energy) ======
  {
    // 7-d input features:{ INTT 3/4 layer R, INTT 3/4 layer Z, INTT 5/4 layer R, INTT 5/4 layer Z, Calo R, Calo Z, Calo Energy }
    std::vector<float> featsMLEproj = { 
      10.0,  5.0,   // INTT 3/4 layer
      15.0,  7.5,   // INTT 5/4 layer
      100.0, 50.0, 8.0 }; // Calo R,Z,Energy
  
    SiCaloPt::InputMLEproj in{featsMLEproj};
    auto r = pTCalc.ComputePt(SiCaloPt::Method::MethodMLEproj, SiCaloPt::AnyInput{in});
    std::cout << "[MLEproj-7D] ok=" << r.ok
              << "  pt=" << r.pt_reco
              << "  err=\"" << r.err << "\"\n";
  }

  // ============ ML：Combined/Gate (2-d input: pt_from_MLEMD, pt_from_MLEproj) ===================
  {
    // 2-d input features:{ pt_from_MLEMD, pt_from_MLEproj }
    std::vector<float> featsMLCombined = {8.0, 9.5};
  
    SiCaloPt::InputMLCombined in{featsMLCombined};
    auto r = pTCalc.ComputePt(SiCaloPt::Method::MethodMLCombined, SiCaloPt::AnyInput{in});
    std::cout << "[MLCombined] ok=" << r.ok
              << "  pt=" << r.pt_reco
              << "  err=\"" << r.err << "\"\n";
  }





  return 0;
};


void analyze_SiSeedPair_pT(
     //const std::string& filename="/sphenix/user/jaein213/tracking/SiliconSeeding/MC/macro/ana/jpsi/ana_all.root", 
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana_jpsi/merged_200k.root", 
//     //const std::string& filename="/sphenix/user/jaein213/tracking/SiliconSeeding/MC/macro/ana/e+/inner_ana_0601_all.root", 
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana_e-/merged_10k.root", 
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana/ana_0_1kevt.root", 
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/mahiro/SIliconCalo/MC/ana/merged_100k.root", // pythia
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_em/all_em.root", // pythia
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_em/all_em1.root", // pythia
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_ep/all_ep.root", // pythia
     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_jpsi/all_jpsi.root", // pythia
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_jpsi/all_jpsi1.root", // pythia
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_jpsi/ana_jpsi_0.root", // pythia
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_em/ana_em_0.root", // pythia
//     const std::string& filename="/sphenix/user/hachiya/INTT/INTT/general_codes/hachiya/SiliconSeeding/macro/data/test_ana2k.root", // pythia
     //const std::string& filename="/sphenix/user/hachiya/INTT/INTT/general_codes/hachiya/SiliconSeeding/macro/data/test_ana5k.root", // pythia
     //const std::string& filename="/sphenix/user/hachiya/INTT/INTT/general_codes/hachiya/SiliconSeeding/macro/data/test_ana10k.root", // pythia
     //const std::string& filename="/sphenix/user/hachiya/INTT/INTT/general_codes/hachiya/SiliconSeeding/macro/data/test_ana50k.root", // pythia
     //const std::string& filename="/sphenix/user/hachiya/INTT/INTT/general_codes/hachiya/SiliconSeeding/macro/data/test_ana.root", // pythia
     //const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana_pythia/ana_addall.root",
//     const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana/all_pythia.root",
//     const std::string& filename="/gpfs/mnt/gpfs02/sphenix/user/hachiya/INTT/INTT/general_codes/hachiya/SiliconSeeding/macro/chargerecalc/ana_test_53879_10evt.root",
     //const std::string& filename="/sphenix/tg/tg01/commissioning/INTT/work/hachiya/SiCaloTrack/data/ana/ana10_pythia.root",
     float phi_threshold = 0.05
)
{

  /////////////////////
  // for PtCalculator
  DemoPaths WS_Path;
  WS_Path.print();
  
  //SiCaloPt::PtCalculatorConfig setting(if you want to use the ML models, Formula-method donnot need these)
  SiCaloPt::PtCalculatorConfig cfg;
  cfg.mlEMD_model_path        = make_opt<decltype(cfg.mlEMD_model_path)      >(WS_Path.emd_onnx);
  cfg.mlEMD_scaler_json       = make_opt<decltype(cfg.mlEMD_scaler_json)     >(WS_Path.emd_scaler_json);
  cfg.mlEproj_model_path      = make_opt<decltype(cfg.mlEproj_model_path)    >(WS_Path.eproj_onnx);
  cfg.mlEproj_scaler_json     = make_opt<decltype(cfg.mlEproj_scaler_json)   >(WS_Path.eproj_scaler_json);
  cfg.mlCombined_model_path   = make_opt<decltype(cfg.mlCombined_model_path) >(WS_Path.combined_onnx);
  cfg.mlCombined_scaler_json  = make_opt<decltype(cfg.mlCombined_scaler_json)>(WS_Path.combined_scaler_json);

  // PtCalculator instance
  //SiCaloPt::PtCalculator calcTutorial(cfg);
  pTCalc.setConfig(cfg);

  // initialize (load models and scalers for ML methods)
  std::string err;
  if (!pTCalc.init(&err)) 
  {
    std::cout << "[init] failed: " << err << std::endl;
    return;
  }
  std::cout << "[init] OK\n";


  demoPt();
  //calcPt();
  //return;
  /////////////////////









  TFile *f = TFile::Open(filename.c_str(), "READ");

  TTree *trackTree = (TTree *)f->Get("trackTree");
  TTree *caloTree  = (TTree *)f->Get("caloTree");
  TTree *siClusTree  = (TTree *)f->Get("SiClusTree");
  TTree *truthTree = (TTree *)f->Get("truthTree");

  int evt;
  std::vector<float> *track_phi    = 0, *track_pt = 0, *track_pz = 0, *track_px = 0, *track_py = 0, *track_eta = 0, *track_z = 0;
  std::vector<float> *track_chi2ndf= 0;
  std::vector<int>   *track_charge = 0, *track_nmaps = 0, *track_nintt = 0;
  std::vector<float> *track_phi_emc= 0, *track_x_emc = 0, *track_y_emc = 0, *track_z_emc = 0;
  std::vector<float> *track_x_oemc = 0, *track_y_oemc = 0, *track_z_oemc = 0;
  std::vector<float> *track_rv_x_emc = 0, *track_rv_y_emc = 0, *track_rv_z_emc = 0;

  int calo_evt;
  std::vector<float> *calo_phi = 0, *calo_energy = 0;
  std::vector<float> *calo_x = 0, *calo_y = 0, *calo_z = 0;
  std::vector<float> *calo_chi2 = 0, *calo_prob = 0;


  std::vector<float>* truth_pt=0, *truth_eta=0, *truth_phi=0;
  std::vector<float> *truth_px = 0, *truth_py = 0, *truth_pz = 0;
  std::vector<int>   *truth_pid=0, *truth_id=0;
  //std::vector<float>* truth_px, *truth_py, *truth_pz, *truth_e;

  // SiClus
  std::vector<int> *SiClus_trackid=0;
  std::vector<int> *SiClus_layer=0;
  std::vector<float> *SiClus_x=0;
  std::vector<float> *SiClus_y=0;
  std::vector<float> *SiClus_z=0;



  trackTree->SetBranchAddress("evt",     &evt);
  trackTree->SetBranchAddress("chi2ndf", &track_chi2ndf);
  trackTree->SetBranchAddress("charge",  &track_charge);
  trackTree->SetBranchAddress("nmaps",   &track_nmaps);
  trackTree->SetBranchAddress("nintt",   &track_nintt);
  trackTree->SetBranchAddress("pt0",     &track_pt);
  trackTree->SetBranchAddress("px0",     &track_px);
  trackTree->SetBranchAddress("py0",     &track_py);
  trackTree->SetBranchAddress("pz0",     &track_pz);
  trackTree->SetBranchAddress("eta0",    &track_eta);
  trackTree->SetBranchAddress("phi0",    &track_phi);
  trackTree->SetBranchAddress("z0",      &track_z);
  trackTree->SetBranchAddress("phi_proj_emc", &track_phi_emc);
  trackTree->SetBranchAddress("x_proj_emc",   &track_x_emc);
  trackTree->SetBranchAddress("y_proj_emc",   &track_y_emc);
  trackTree->SetBranchAddress("z_proj_emc",   &track_z_emc);
  trackTree->SetBranchAddress("x_proj_oemc",  &track_x_oemc);
  trackTree->SetBranchAddress("y_proj_oemc",  &track_y_oemc);
  trackTree->SetBranchAddress("z_proj_oemc",  &track_z_oemc);
  trackTree->SetBranchAddress("x_rv_proj_emc",  &track_rv_x_emc);
  trackTree->SetBranchAddress("y_rv_proj_emc",  &track_rv_y_emc);
  trackTree->SetBranchAddress("z_rv_proj_emc",  &track_rv_z_emc);
/*
  trackTree->SetBranchAddress("pt0",      &track_pt);
  trackTree->SetBranchAddress("eta0",     &track_eta);
  trackTree->SetBranchAddress("phi0",     &track_phi);
  trackTree->SetBranchAddress("z0",       &track_z);
  trackTree->SetBranchAddress("phi_proj_emc", &track_phi_emc);
  trackTree->SetBranchAddress("z_proj_emc",   &track_z_emc);
*/

///////////
   // Truth tree and branches
  if(truthTree!=nullptr){
    truthTree->SetBranchAddress("truth_pt",  &truth_pt);
    truthTree->SetBranchAddress("truth_pz",  &truth_pz);
    truthTree->SetBranchAddress("truth_eta", &truth_eta);
    truthTree->SetBranchAddress("truth_phi", &truth_phi);
    truthTree->SetBranchAddress("truth_px",  &truth_px);
    truthTree->SetBranchAddress("truth_py",  &truth_py);
    truthTree->SetBranchAddress("truth_pid", &truth_pid);
    truthTree->SetBranchAddress("truth_id",  &truth_id);
  }

  //--truthTree->SetBranchAddress("truth_pid", &truth_pid);
  //--truthTree->SetBranchAddress("truth_px",  &truth_px);
  //--truthTree->SetBranchAddress("truth_py",  &truth_py);
  //--truthTree->SetBranchAddress("truth_e",   &truth_e);
  //--truthTree->SetBranchAddress("truth_pt",  &truth_pt);

///////////
  siClusTree->SetBranchAddress("Siclus_trackid", &SiClus_trackid);
  siClusTree->SetBranchAddress("Siclus_layer",   &SiClus_layer);
  siClusTree->SetBranchAddress("Siclus_x",       &SiClus_x);
  siClusTree->SetBranchAddress("Siclus_y",       &SiClus_y);
  siClusTree->SetBranchAddress("Siclus_z",       &SiClus_z);


///////////
  caloTree->SetBranchAddress("calo_evt", &calo_evt);
  caloTree->SetBranchAddress("phi",      &calo_phi);
  caloTree->SetBranchAddress("x",        &calo_x);
  caloTree->SetBranchAddress("y",        &calo_y);
  caloTree->SetBranchAddress("z",        &calo_z);
  caloTree->SetBranchAddress("energy",   &calo_energy);
  caloTree->SetBranchAddress("chi2",     &calo_chi2);
  caloTree->SetBranchAddress("prob",     &calo_prob);

///////////

  //TFile *outFile = new TFile("dphi_distribution_em.root", "RECREATE");
  TFile *outFile = new TFile("dphi_distribution.root", "RECREATE");
  TH1F *h_dphi = new TH1F("h_dphi", "Track - Calo #Delta#phi;#Delta#phi;Counts", 200, -0.3, 0.3);
  TH1F *h_dphi_emc = new TH1F("h_dphi (EMCal Proj)", "Track - Calo #Delta#phi;#Delta#phi;Counts", 200, -0.3, 0.3);
  TH2F *h_dphi_emc_pt = new TH2F("h_dphi_pt", "Track - Calo #Delta#phi;pT;#Delta#phi", 200, 0, 20, 200, -0.3, 0.3);
  TH2F *h_dphi_emc_pt_truth = new TH2F("h_dphi_pt_truth", "Track - Calo #Delta#phi;pT{true};#Delta#phi", 200, 0, 20, 200, -0.3, 0.3);
  TH1F *h_EoverP_all = new TH1F("h_EoverP_all", "E/p of matched track-calo;E/p;Counts", 100, 0, 5);
  TH1F *h_EoverP_cut = new TH1F("h_EoverP_cut", "E/p of matched track-calo phi cut;E/p;Counts", 100, 0, 5);
  TH1F *h_dz = new TH1F("h_dz", "Track - Calo #Delta z;#Delta z (cm);Counts", 200, -50, 50);
  TH1F *h_dz_emc = new TH1F("h_dz_emc", "Track(EMCal Proj) - Calo #Delta z;#Delta z (cm);Counts", 200, -50, 50);

  TH1F *h_mass     = new TH1F("h_mass", "Invariant mass of matched track pairs (e^{+}e^{-});Mass (GeV/c^{2});Counts", 200, 0, 5);
  // Invariant mass histograms for different J/psi pT bins
  TH1F *h_mass_pt0_1 = new TH1F("h_mass_pt0_1", "Invariant mass (pT 0-1 GeV/c);Mass (GeV/c^{2});Counts", 200, 0, 5);
  TH1F *h_mass_pt1_2 = new TH1F("h_mass_pt1_2", "Invariant mass (pT 1-2 GeV/c);Mass (GeV/c^{2});Counts", 200, 0, 5);
  TH1F *h_mass_pt2_3 = new TH1F("h_mass_pt2_3", "Invariant mass (pT 2-3 GeV/c);Mass (GeV/c^{2});Counts", 200, 0, 5);
  TH1F *h_mass_pt3_4 = new TH1F("h_mass_pt3_4", "Invariant mass (pT 3-4 GeV/c);Mass (GeV/c^{2});Counts", 200, 0, 5);
  TH1F *h_mass_pt4up = new TH1F("h_mass_pt4up", "Invariant mass (pT >4 GeV/c);Mass (GeV/c^{2});Counts", 200, 0, 5);
  TH1F *h_track_chi2ndf_matched = new TH1F("h_track_chi2ndf_matched", "Chi2/NDF of matched tracks;#chi^{2}/ndf;Counts", 100, 0, 10);

  TH2F *h_reco_vs_truth_pt = new TH2F("h_reco_vs_truth_pt", "Reconstructed pT vs Truth pT;Truth pT (GeV/c);Reconstructed pT (GeV/c)", 100, 0, 20, 100, 0, 20);

  TNtuple *h_ntp_sicalo = new TNtuple("ntp_sicalo", "SiSeed + Calo combination",       "pt:pz:c:chi2:nintt:nmaps:hitbit:eta:the:phi:z:phi_pemc:z_pemc:phi_poemc:z_poemc:phi_rvpemc:z_rvpemc:phi_intt:pt_tr:pz_tr:phi_tr:pid_tr:nemc:phi_emc:z_emc:e_emc:chi2_emc:phi_calo:pt_calo");

  TNtuple *h_ntp_sicalotrk = new TNtuple("ntp_sicalotrk", "SiSeed + Calo combination", "pt:pz:c:chi2:nintt:nmaps:hitbit:eta:the:phi:z:phi_pemc:z_pemc:phi_poemc:z_poemc:phi_rvpemc:z_rvpemc:phi_intt:pt_tr:pz_tr:phi_tr:pid_tr:nemc:phi_emc:z_emc:e_emc:chi2_emc:phi_calo:pt_calo");

  TNtuple *h_ntp_pair = new TNtuple("ntp_pair", "pair", "mass:pt:phi:eta:massc:ptc:ptp:pzp:eopp:dpp:dzp:ptm:pzm:eopm:dpm:dzm");


  // Define a struct to hold matched track info
  struct MatchedTrack {
    size_t track_idx;
    size_t calo_idx;
    float  min_dphi;
    float  min_dphi_emc;
    float  min_dz;
    float  min_dz_emc;
    float  eop;
    float  pt_calo;
    float  eop_calo;
    int    nmaps;
    int    nintt;
    int    charge;
    float  chi2ndf;
  };

  struct stCluster {
    int   m_clsTrkid;
    int   m_clslyr  ;
    float m_clsx    ;
    float m_clsy    ;
    float m_clsz    ;
  };
  ////////////////////////////////////////////////////

  Long64_t nentries = trackTree->GetEntries();

  int nskip=0;

  for (Long64_t i = 0; i < nentries; ++i) // Event Loop
  {

    trackTree->GetEntry(i);
    caloTree->GetEntry(i);
    siClusTree->GetEntry(i);
    if(truthTree!=nullptr) truthTree->GetEntry(i);

    if (evt != calo_evt)
    {
      std::cerr << "Warning: evt mismatch at entry " << i << ": track evt = " << evt << ", calo evt = " << calo_evt << std::endl;
      continue;
    }

    if((evt%1000)==0) {
       std::cout << "Matching evt = " << evt << std::endl;
    }

    std::vector<stCluster> vClus;



    int nclusSize = SiClus_trackid->size();

    cout<<evt <<" : ntrk, nemc : "<<track_phi->size()<<", "<< calo_phi->size()<< " "<<nclusSize<<endl;
    std::vector<MatchedTrack> matched_tracks;

    int itotalClus=0;
    ///////////////////////////
    // First pass: find matched tracks and store info
    for (size_t it = 0; it < track_phi->size(); ++it)
    {
      int nmaps = (*track_nmaps)[it];
      int nintt = (*track_nintt)[it];

      float& pt     = ((*track_pt)[it]);
      float& px     = ((*track_px)[it]);
      float& py     = ((*track_py)[it]);
      float& pz     = ((*track_pz)[it]);
      float& eta    = ((*track_eta)[it]);
      int&   charge = ((*track_charge)[it]);

      int nclus = nmaps+nintt;
      std::cout << "charge : " << charge<<" "<<pt<<" "<<nmaps<<" "<<nintt<<" "<<calo_phi->size()<<" "<<nclus<<std::endl;

      // cluster:
      int hitbit=0;
      stCluster iCl, oCl;
      for(int ic=0; ic<nclus; ic++) {
        int   clsid    = ic+itotalClus;

        int layer =  ((*SiClus_layer)[clsid]);
        stCluster Cl {
          ((*SiClus_trackid)[clsid]),
          layer,
          ((*SiClus_x)[clsid]),
          ((*SiClus_y)[clsid]),
          ((*SiClus_z)[clsid]) 
        };
        vClus.push_back(Cl);

        if(0<=Cl.m_clslyr&&Cl.m_clslyr<5) iCl = Cl;
        else                              oCl = Cl;

        hitbit |= (1<<layer);

        //cout<<"cls : "<<clsid<<" "<<Cl.m_clsTrkid<<" "<<Cl.m_clslyr<<" "<<Cl.m_clsx<<" "<<Cl.m_clsy<<" "<<Cl.m_clsz<<endl;
      }
      itotalClus+=nclus;

      // truth track
      float pt_tr=0, pz_tr=0, phi_tr=0, eta_tr=0;
      int pid_tr=0;

      if(truthTree!=nullptr){
        int ntruth = truth_pt->size();
        float min_dp = 99999; 
        size_t min_itr=1000000;
        for (size_t itr = 0; itr < ntruth; ++itr)
        {
          int pid = ((*truth_pid)[itr]);
          //if( !(abs(pid)==11 // e-
          //  || abs(pid)==13 // mu-
          //  || abs(pid)==211 // pi+
          //  || abs(pid)==321 // K+
          //  || abs(pid)==2212 // p
          //  )) continue;


          float dpx = (px - ((*truth_px)[itr]));
          float dpy = (py - ((*truth_py)[itr]));
          float dpz = (pz - ((*truth_pz)[itr]));
          float dp  = sqrt(dpx*dpx+dpy*dpy+dpz*dpz);
          if(dp<min_dp) {
            min_dp = dp;
            min_itr = itr;
          }

          //--cout<<"tru itr: "<<itr<<" "<<pid<<" "<<dp<<endl;
          //--pt_tr  = ((*truth_pt)[itr]);
          //--pz_tr  = ((*truth_pz)[itr]);
          //--phi_tr = ((*truth_phi)[itr]);
          //--eta_tr = ((*truth_eta)[itr]);
          //--px_tr  = ((*truth_px)[itr]);
          //--py_tr  = ((*truth_py)[itr]);
        }

        if(min_itr<ntruth){
          pt_tr  = ((*truth_pt)[min_itr]);
          pz_tr  = ((*truth_pz)[min_itr]);
          phi_tr = ((*truth_phi)[min_itr]);
          eta_tr = ((*truth_eta)[min_itr]);
          pid_tr = ((*truth_pid)[min_itr]);
        }
        cout<<"tru : "<<ntruth<<" "<<min_itr<<" "<<pt_tr<<" "<<pz_tr<<" "<<phi_tr<<" "<<eta_tr<<" "<<pid_tr<<endl;
      }

      if (pt < 0.3) continue;

      if (nmaps < 1 || nintt < 1)
      {
        if(nskip<1000) {
          std::cout << "Skipping track with nmaps = " << nmaps << " and nintt = " << nintt << std::endl;
        } else if(nskip==1000) {
          std::cout << "exceed nskip. comment suppressed" << std::endl;
        }
        nskip++;

        continue;
      }

      //cout<<"cls-in : "<<iCl.m_clsTrkid<<" "<<iCl.m_clslyr<<" "<<iCl.m_clsx<<" "<<iCl.m_clsy<<" "<<iCl.m_clsz<<endl;
      //cout<<"cls-out: "<<oCl.m_clsTrkid<<" "<<oCl.m_clslyr<<" "<<oCl.m_clsx<<" "<<oCl.m_clsy<<" "<<oCl.m_clsz<<endl;
      //std::cout << "charge : " << charge<<" "<<pt<<std::endl;

      float phi_intt = atan2(oCl.m_clsy-iCl.m_clsy, oCl.m_clsx-iCl.m_clsx);

      float& x_emc = (*track_x_emc)[it];   //12
      float& y_emc = (*track_y_emc)[it];   //12
      float phi_proj = atan2(y_emc, x_emc);

      float& x_oemc = (*track_x_oemc)[it];   //12
      float& y_oemc = (*track_y_oemc)[it];   //12
      float phi_proj_o = atan2(y_oemc, x_oemc);

      float& x_rv_emc = (*track_rv_x_emc)[it];   //12
      float& y_rv_emc = (*track_rv_y_emc)[it];   //12
      float phi_proj_rv = atan2(y_rv_emc, x_rv_emc);

      // Find calo cluster with minimum dphi_emc for this track
      float ntp_value[30] = {
                pt,                   //0
                pz,                   //1
                (float)(charge),      //2
                (*track_chi2ndf)[it], //3
                (float)nintt,         //4
                (float)nmaps,         //5
                (float)hitbit,        //6
                (*track_eta)[it],     //7
                (float)(2.* std::atan( std::exp(-eta) )),//8
                (*track_phi)[it],     //9   
                (*track_z)[it],       //10
                phi_proj, //(*track_phi_emc)[it], //11
                (*track_z_emc)[it],   //12
                phi_proj_o, //(*track_phi_emc)[it], //13
                (*track_z_oemc)[it],   //14
                phi_proj_rv, //(*track_phi_emc)[it], //15
                (*track_rv_z_emc)[it],   //16
                phi_intt,             //17
                pt_tr,                //18
                pz_tr,                //19
                phi_tr,               //20
                (float)pid_tr,        //21
                (float)calo_phi->size() //22
              };

      PtInput input;
      input.SiIn[0] = iCl.m_clsx;
      input.SiIn[1] = iCl.m_clsy;
      input.SiIn[2] = iCl.m_clsz;
      input.SiOut[0] = oCl.m_clsx;
      input.SiOut[1] = oCl.m_clsy;
      input.SiOut[2] = oCl.m_clsz;

      //const float par[2] = {0.21, -0.986}; // par for pT_calo calculation
      const float par[2] = {0.2, -1}; // par for pT_calo calculation

      float min_dr = 1e9;
      size_t min_ic = calo_phi->size();
      float min_dphi = 0, min_dphi_emc = 0, min_dz = 0, min_dz_emc = 0;
      float match_pt_calo=0;
      //cout<<"aaaa"<<endl;
      for (size_t ic = 0; ic < calo_phi->size(); ++ic)
      {
      //--cout<<"bbb"<<endl;
       
        //float dphi_emc = (*track_phi_emc)[it] - (*calo_phi)[ic];
        float dphi_emc = (*calo_phi)[ic] - phi_proj;
        if (dphi_emc >  M_PI) dphi_emc -= 2 * M_PI;
        if (dphi_emc < -M_PI) dphi_emc += 2 * M_PI;

        const float p0 = -0.181669;
        const float p1 =  0.00389827;
        //float dphi_emc_corr = dphi_emc - charge*(p0/pt + p1);
        float dphi_emc_corr = charge*dphi_emc - 0.18*std::pow(pt, -0.986);

        float x_calo = (*calo_x)[ic];
        float y_calo = (*calo_y)[ic];
        float z_calo = (*calo_z)[ic];

        float dz_emc = z_calo - (*track_z_emc)[it];


        float phi_calo = atan2(y_calo - oCl.m_clsy,  x_calo - oCl.m_clsx);

        float dphi = phi_calo - phi_intt;
        //float pt_calo = par[0]*pow(-charge*dphi, par[1]);//cal_CaloPt(dphi);

        float pt_inv = charge*0.02+4.9*(-charge*dphi)-0.6*pow(-charge*dphi, 2);
        float pt_calo = 1./pt_inv;

        input.Calo[0] = x_calo;
        input.Calo[1] = y_calo;
        input.Calo[2] = z_calo;
        input.CaloE   = (*calo_energy)[ic];

         
        //pt_calo = calcPt(input);

        ntp_value[23] = (*calo_phi)[ic];
        ntp_value[24] = z_calo;
        ntp_value[25] = (*calo_energy)[ic];
        ntp_value[26] = (*calo_chi2)[ic];
        ntp_value[27] = phi_calo;
        ntp_value[28] = pt_calo;
        
        h_ntp_sicalo->Fill(ntp_value);
        //--std::cout << "   energy : " << ntp_value[10]<<std::endl;

        float dphi_emc_cm = dphi_emc*93.5/3.; // normalized by sigma=3
        float dr = sqrt(dphi_emc_cm*dphi_emc_cm + dz_emc*dz_emc);

        //if (fabs(dphi_emc) < min_dr)
        //if (fabs(dphi_emc_corr) < min_dr)
        //if (fabs(dz_emc) < min_dr)
        if (fabs(dr) < min_dr)
        {
          //min_dr = fabs(dphi_emc);
          //min_dr = fabs(dphi_emc_corr);
          
          //min_dr = fabs(dz_emc);
          min_dr = fabs(dr);
          min_ic = ic;

          //min_dphi = (*track_phi)[it] - (*calo_phi)[ic];
          //if (min_dphi >  M_PI) min_dphi -= 2 * M_PI;
          //if (min_dphi < -M_PI) min_dphi += 2 * M_PI;

          //min_dphi_emc = dphi_emc_corr;
          //min_dz = (*track_z)[it] - (*calo_z)[ic];
          //min_dz_emc = (*track_z_emc)[it] - (*calo_z)[ic];
          min_dphi_emc = dphi_emc;
          min_dz_emc   = dz_emc;

          match_pt_calo = pt_calo;
        }
      }

      if (min_ic != calo_phi->size())
      { // for matched track
        float x_calo = (*calo_x)[min_ic];
        float y_calo = (*calo_y)[min_ic];
        float z_calo = (*calo_z)[min_ic];

        float phi_calo = atan2(y_calo - oCl.m_clsy,  x_calo - oCl.m_clsx);

        float dphi = phi_calo - phi_intt;
        //float pt_calo = par[0]*pow(-charge*dphi, par[1]);//cal_CaloPt(dphi);
        float pt_inv = charge*0.02+4.9*(-charge*dphi)-0.6*pow(-charge*dphi, 2);
        float pt_calo = 1./pt_inv;

        cout<<"new pT : "<<pt_calo<<endl;

        input.Calo[0] = x_calo;
        input.Calo[1] = y_calo;
        input.Calo[2] = z_calo;
        input.CaloE   = (*calo_energy)[min_ic];
        calcPt(input);
        //pt_calo = calcPt(input);

        ntp_value[23] = (*calo_phi)[min_ic];
        ntp_value[24] = z_calo;
        ntp_value[25] = (*calo_energy)[min_ic];
        ntp_value[26] = (*calo_chi2)[min_ic];
        ntp_value[27] = phi_calo;
        ntp_value[28] = pt_calo;
      }
      else {
        ntp_value[23] = -9999.;
        ntp_value[24] = -9999.;
        ntp_value[25] = -9999.;
        ntp_value[26] = -9999.;
        ntp_value[27] = -9999.;
        ntp_value[28] = -9999.;
      }
      h_ntp_sicalotrk->Fill(ntp_value);

      if (min_ic == calo_phi->size()) {
        continue; // no match
      }

      float p = (*track_pt)[it] * std::cosh((*track_eta)[it]);
      float E = (*calo_energy)[min_ic];
      float eop = -999;
      if (p > 0)
      {
        eop = E / p;
        h_EoverP_all->Fill(eop);
      }
      float p_calo = match_pt_calo * std::cosh((*track_eta)[it]);
      float eop_calo = E/p_calo;

      h_dphi->Fill(min_dphi);
      h_dphi_emc->Fill(min_dphi_emc);
      h_dphi_emc_pt->Fill((*track_pt)[it], min_dphi_emc);
      //h_dphi_emc_pt_truth->Fill((*truth_pt)[0], min_dphi_emc);
      h_dz->Fill(min_dz);
      h_dz_emc->Fill(min_dz_emc);
     // std::cout << "  ➤ Δφ = " << min_dphi << ", Δz = " << min_dz << std::endl;

      h_track_chi2ndf_matched->Fill((*track_chi2ndf)[it]);
      //--if ((*track_chi2ndf)[it] > 5)
      //--  continue; // Skip tracks with high chi2/ndf
      if (min_dz_emc > 4)
        continue;
     // h_reco_vs_truth_pt->Fill((*truth_pt)[0], (*track_pt)[it]);

      matched_tracks.push_back({it, min_ic, min_dphi, min_dphi_emc, min_dz, min_dz_emc, eop, match_pt_calo, eop_calo, nmaps, nintt, (*track_charge)[it], (*track_chi2ndf)[it]});
    }

    //cout<<"nclus-end : "<<itotalClus<<endl;

    /////////////////////////
    //--// 1.5 pass: split positive and negative to different array;
    //--std::vector<MatchedTrack> v_trkp, v_trkm;
    //--for (size_t idx = 0; idx < matched_tracks.size(); ++idx)
    //--{
    //--  const auto& mt = matched_tracks[idx];
    //--  if(mt.charge>0) v_trkp.push_back(mt);
    //--  else            v_trkm.push_back(mt);
    //--}

    /////////////////////////
    // Second pass: loop through matched tracks and fill histograms and compute invariant mass
    cout<<"matched track : "<<matched_tracks.size()<<endl;
    float v_pair[20];
    for (size_t idx = 0; idx < matched_tracks.size(); ++idx)
    {
      const auto& mt = matched_tracks[idx];

      if (mt.nmaps > 2 && mt.nintt > 1 && std::abs(mt.min_dz_emc) < 4)
      {
        h_EoverP_cut->Fill(mt.eop);

        //--if (!(mt.eop > 0.01 && mt.eop < 2)) continue;

        for (size_t jdx = idx+1; jdx < matched_tracks.size(); ++jdx)
        {
          const auto& mt2 = matched_tracks[jdx];

          //if (mt2.nmaps < 2 || mt2.nintt < 1) continue;
          //if (std::abs(mt2.min_dz_emc) > 4) continue;

          if (mt2.nmaps > 2 && mt2.nintt > 1 && std::abs(mt2.min_dz_emc) < 4) {

            if (mt.charge * mt2.charge >= 0) continue; // opposite sign only

            //float eop2 = -999;
            //float p2 = (*track_pt)[mt2.track_idx] * std::cosh((*track_eta)[mt2.track_idx]);
            //float E2_raw = (*calo_energy)[mt2.calo_idx];
            //if (p2 > 0) eop2 = E2_raw / p2;
            //if (!(eop2 > 0.8 && eop2 < 2)) continue;

            //--if (!(mt2.eop > 0.6 && mt2.eop < 2)) continue;

            float pt1 = (*track_pt)[mt.track_idx],  eta1 = (*track_eta)[mt.track_idx],  phi1 = (*track_phi)[mt.track_idx];
            float pt2 = (*track_pt)[mt2.track_idx], eta2 = (*track_eta)[mt2.track_idx], phi2 = (*track_phi)[mt2.track_idx];
            TLorentzVector v1, v2;
            v1.SetPtEtaPhiM(pt1, eta1, phi1, electron_mass);
            v2.SetPtEtaPhiM(pt2, eta2, phi2, electron_mass);
            //v1.SetPtEtaPhiM(pt1, eta1, phi1, pion_mass);
            //v2.SetPtEtaPhiM(pt2, eta2, phi2, pion_mass);

            TLorentzVector vpair = v1 + v2;
 
            //float mass = (v1 + v2).M();
            //float total_pt = pt1 + pt2;
            float mass     = vpair.M();
            float total_pt = vpair.Pt();
        
            // mass with pt_calo
            float ptc1 = mt.pt_calo;
            float ptc2 = mt2.pt_calo;
            TLorentzVector vc1, vc2;
            vc1.SetPtEtaPhiM(ptc1, eta1, phi1, electron_mass);
            vc2.SetPtEtaPhiM(ptc2, eta2, phi2, electron_mass);
            //vc1.SetPtEtaPhiM(ptc1, eta1, phi1, pion_mass);
            //vc2.SetPtEtaPhiM(ptc2, eta2, phi2, pion_mass);
            TLorentzVector vpairc = vc1 + vc2;
 
            float massc     = vpairc.M();
            float total_ptc = vpairc.Pt();

            //if ((0.8 < mt.eop  && mt.eop  < 2) &&
            //    (0.8 < mt2.eop && mt2.eop < 2)) 
            {
              h_mass->Fill(mass);

              // Fill the appropriate mass histogram based on total pT
              if (total_pt < 1.0) {
                h_mass_pt0_1->Fill(mass);
              } else if (total_pt < 2.0) {
                h_mass_pt1_2->Fill(mass);
              } else if (total_pt < 3.0) {
                h_mass_pt2_3->Fill(mass);
              } else if (total_pt < 4.0) {
                h_mass_pt3_4->Fill(mass);
              } else {
                h_mass_pt4up->Fill(mass);
              }
            }

    //TNtuple *h_ntp_pair = new TNtuple("ntp_pair", "SiSeed + Calo combination", 
    //"mass:pt:phi:eta:ptp:pzp:phip:dpp:dzp:ptm:pzm:phim:dpm:dzm");
 
            v_pair[ 0] = mass;
            v_pair[ 1] = total_pt;
            v_pair[ 2] = vpair.Phi();
            v_pair[ 3] = vpair.Eta();
            v_pair[ 4] = massc; // mass with pt_calo
            v_pair[ 5] = total_ptc; // pair-pt with pt_calo
            v_pair[ 6] = pt1;
            v_pair[ 7] = v1.Pz();
            v_pair[ 8] = mt.eop; //mt.eop_calo;
            v_pair[ 9] = mt.min_dphi_emc;
            v_pair[10] = mt.min_dz_emc;
            v_pair[11] = pt2;
            v_pair[12] = v2.Pz();
            v_pair[13] = mt2.eop; //mt2.eop_calo;
            v_pair[14] = mt2.min_dphi_emc;
            v_pair[15] = mt2.min_dz_emc;

            h_ntp_pair->Fill(v_pair);
          }
        }
      }
    }
  }

  std::cout<<"total nskip : "<< nskip <<std::endl;

 // TFile *outFile = new TFile("dphi_distribution_e-.root", "RECREATE");
  //TFile *outFile = new TFile("dphi_distribution_PYTHIA_temp.root", "RECREATE");
  std::cout << "Saving histogram to dphi_distribution.root" << std::endl;
  h_dphi->Write();
  h_dphi_emc->Write();
  h_dphi_emc_pt->Write();
  h_dphi_emc_pt_truth->Write();
  h_EoverP_all->Write();
  h_EoverP_cut->Write();
  h_dz->Write();
  h_dz_emc->Write();
  h_mass->Write();
  h_mass_pt0_1->Write();
  h_mass_pt1_2->Write();
  h_mass_pt2_3->Write();
  h_mass_pt3_4->Write();
  h_mass_pt4up->Write();
  h_track_chi2ndf_matched->Write();
  h_reco_vs_truth_pt->Write();
  h_ntp_sicalo->Write();
  h_ntp_sicalotrk->Write();
  h_ntp_pair->Write();

  outFile->Close();

  std::cout << "Δφ histogram saved to 'dphi_distribution.root'" << std::endl;
  f->Close();
}
