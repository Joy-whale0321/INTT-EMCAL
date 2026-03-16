#include <TSystem.h>
#include <TStopwatch.h>
#include <iostream>
#include <vector>

// #include "PtCalculator.h"  // SiCaloPt::PtCalculator & friends
#include "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/src/PtCalculator.h"  // SiCaloPt::PtCalculator
// #include "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/src/PtCalculator.h"  // SiCaloPt::PtCalculator

R__LOAD_LIBRARY(/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/src/libPtCalc.so)
// R__LOAD_LIBRARY(/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/src/libPtCalc.so)


// ---- Weights(onnx) and Scalers(json) Path on rcf ---------------------------
struct DemoPaths
{
    std::string emd_onnx_m10         = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_eP_s-10.onnx"; 
    std::string emd_scaler_json_m10  = "";

    std::string emd_onnx_m11         = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_eI_s-11.onnx"; 
    std::string emd_scaler_json_m11  = "";

    std::string emd_onnx_m12         = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_eD_s-12.onnx"; 
    std::string emd_scaler_json_m12  = "";

    std::string emd_onnx_p10         = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_pP_s+10.onnx"; 
    std::string emd_scaler_json_p10  = "";

    std::string emd_onnx_p11         = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_pI_s+11.onnx"; 
    std::string emd_scaler_json_p11  = "";

    std::string emd_onnx_p12         = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_pD_s+12.onnx"; 
    std::string emd_scaler_json_p12  = "";

    std::string eproj_onnx        = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEproj.onnx"; 
    std::string eproj_scaler_json = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/scaler_MLEproj.json";

    std::string combined_onnx         = "/sphenix/user/jzhang1/testcode4all/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLCombined.onnx"; 
    std::string combined_scaler_json  = "";
};


// ---- Weights(onnx) and Scalers(json) Path local ---------------------------
// struct DemoPaths
// {
//     std::string emd_onnx_m10         = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_eP_s-10.onnx"; 
//     std::string emd_scaler_json_m10  = "";

//     std::string emd_onnx_m11         = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_eI_s-11.onnx"; 
//     std::string emd_scaler_json_m11  = "";

//     std::string emd_onnx_m12         = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_eD_s-12.onnx"; 
//     std::string emd_scaler_json_m12  = "";

//     std::string emd_onnx_p10         = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_pP_s+10.onnx"; 
//     std::string emd_scaler_json_p10  = "";

//     std::string emd_onnx_p11         = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_pI_s+11.onnx"; 
//     std::string emd_scaler_json_p11  = "";

//     std::string emd_onnx_p12         = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEMD_pD_s+12.onnx"; 
//     std::string emd_scaler_json_p12  = "";

//     std::string eproj_onnx        = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLEproj.onnx"; 
//     std::string eproj_scaler_json = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/scaler_MLEproj.json";

//     std::string combined_onnx         = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ML4Reco/Implement/ML_Weight_Scaler/model_MLCombined.onnx"; 
//     std::string combined_scaler_json  = "";
// };

// ---- turn string into optional<string> for Config ------------------------
template<typename Opt>
Opt make_opt(const std::string& s) 
{
    if (s.empty()) return std::nullopt;
    return s;
}

// ---- PtCalc Tutorial -------------------------------------
void PtCalcMLTutorial()
{
    // Use appropriate paths in your environment, the default "DemoPaths" setup is correct here for mine
    DemoPaths WS_Path;  

    // SiCaloPt::PtCalculatorConfig setting(if you want to use the ML models, Formula-method donnot need these)
    SiCaloPt::PtCalculatorConfig cfg;
    cfg.mlEMD_model_path[-10]        = WS_Path.emd_onnx_m10;
    cfg.mlEMD_model_path[-11]        = WS_Path.emd_onnx_m11;
    cfg.mlEMD_model_path[-12]        = WS_Path.emd_onnx_m12;
    cfg.mlEMD_model_path[10]         = WS_Path.emd_onnx_p10;
    cfg.mlEMD_model_path[11]         = WS_Path.emd_onnx_p11;
    cfg.mlEMD_model_path[12]         = WS_Path.emd_onnx_p12;

    cfg.mlEproj_model_path      = make_opt<decltype(cfg.mlEproj_model_path)>(WS_Path.eproj_onnx);
    cfg.mlEproj_scaler_json     = make_opt<decltype(cfg.mlEproj_scaler_json)>(WS_Path.eproj_scaler_json);
    cfg.mlCombined_model_path   = make_opt<decltype(cfg.mlCombined_model_path)>(WS_Path.combined_onnx);
    cfg.mlCombined_scaler_json  = make_opt<decltype(cfg.mlCombined_scaler_json)>(WS_Path.combined_scaler_json);

    // PtCalculator instance
    SiCaloPt::PtCalculator calcTutorial(cfg);

    // set cluster reco mode, 0 = Projected(fixed radius 93.5cm), 1 = InnerFace center, 2 = geomtric center
    calcTutorial.setClusterRecoMode(0); 
    
    // initialize (load models and scalers for ML methods)
    std::string err;
    if (!calcTutorial.init(&err)) 
    {
        std::cout << "[init] failed: " << err << std::endl;
        return;
    }
    std::cout << "[init] OK\n";
 
    // auto r is PtResult struct, with r.pt_reco is the reconstructed pt, r.ok is whether successful, and r.err is error message if any.
    // ======================== EMD Formula ==========================
    {
        SiCaloPt::InputEMD in;
        in.EMD_Angle  = -0.20;  // delta_phi - EM Deflection angle in rad
        in.EMD_Eta    = 0.5;   // EMCal Cluster eta

        auto r = calcTutorial.ComputePt(SiCaloPt::Method::MethodEMD, SiCaloPt::AnyInput{in});
        std::cout << "[EMD-analytic] ok=" << r.ok
                  << "  pt=" << r.pt_reco
                  << "  err=\"" << r.err << "\"\n";
    }

    // ======================== Eproj Formula ========================
    {
        SiCaloPt::InputEproj in;
        in.Energy_Calo   = 1.8;   // EMCal Cluster energy in GeV
        in.Radius_Calo   = 93.5;  // EMCal Cluster radius in cm
        in.Z_Calo        = 0.0;   // EMCal Cluster z in cm
        in.Radius_vertex = 0.0;   // Vertex radius in cm
        in.Z_vertex      = 0.0;   // Vertex z in cm

        auto r = calcTutorial.ComputePt(SiCaloPt::Method::MethodEproj, SiCaloPt::AnyInput{in});
        std::cout << "[Eproj-analytic] ok=" << r.ok
                  << "  pt=" << r.pt_reco
                  << "  err=\"" << r.err << "\"\n";
    }

    // ============= ML：MLEMD (2-d input: dphi_EMD, eta_track) ===============
    {
        // 2-d input features:{ dphi_EMD, eta_track }
        std::vector<float> featsMLEMD = {-0.20, 0.5};

        SiCaloPt::InputMLEMD in{featsMLEMD};
        auto r = calcTutorial.ComputePt(SiCaloPt::Method::MethodMLEMD, SiCaloPt::AnyInput{in});
        std::cout << "[MLEMD-2D] ok=" << r.ok
                << "  pt=" << r.pt_reco
                << "  err=\"" << r.err << "\"\n";
    }

    // ====== ML：MLEproj (7-d input: INTT 3/4 layer R,Z, INTT 5/4 layer R,Z, Calo R,Z,Energy) ======
    {
        // 7-d input features:{ INTT 3/4 layer R, INTT 3/4 layer Z, INTT 5/4 layer R, INTT 5/4 layer Z, Calo R, Calo Z, Calo Energy }
        std::vector<float> featsMLEproj = { 10.0,  5.0,   // INTT 3/4 layer
                                            15.0,  7.5,   // INTT 5/4 layer
                                            100.0, 50.0, 8.0 }; // Calo R,Z,Energy

        SiCaloPt::InputMLEproj in{featsMLEproj};
        auto r = calcTutorial.ComputePt(SiCaloPt::Method::MethodMLEproj, SiCaloPt::AnyInput{in});
        std::cout << "[MLEproj-7D] ok=" << r.ok
                << "  pt=" << r.pt_reco
                << "  err=\"" << r.err << "\"\n";
    }

    // ============ ML：Combined/Gate (2-d input: pt_from_MLEMD, pt_from_MLEproj) ===================
    {
        // 2-d input features:{ pt_from_MLEMD, pt_from_MLEproj }
        std::vector<float> featsMLCombined = {8.0, 9.5};

        SiCaloPt::InputMLCombined in{featsMLCombined};
        auto r = calcTutorial.ComputePt(SiCaloPt::Method::MethodMLCombined, SiCaloPt::AnyInput{in});
        std::cout << "[MLCombined] ok=" << r.ok
                << "  pt=" << r.pt_reco
                << "  err=\"" << r.err << "\"\n";
    }

    // // ======================== Batch example ========================
    // {
    //     std::vector<SiCaloPt::InputEproj> batch(1000);
    //     for (size_t i = 0; i < batch.size(); ++i) {
    //     auto& x = batch[i];
    //     x.Energy_Calo   = 1.0 + 0.001 * i;
    //     x.Radius_Calo   = 93.5;
    //     x.Z_Calo        = 0.0;
    //     x.Radius_vertex = 0.0;
    //     x.Z_vertex      = 0.0;
    //     }

    //     TStopwatch sw; sw.Start();
    //     double sum_pt = 0.0; size_t ok_cnt = 0;
    //     for (auto& x : batch) {
    //     auto r = calcTutorial.ComputePt(SiCaloPt::Method::MethodEproj, SiCaloPt::AnyInput{x});
    //     if (r.ok) { sum_pt += r.pt_reco; ++ok_cnt; }
    //     }
    //     sw.Stop();

    //     std::cout << "[Batch Eproj] ok=" << ok_cnt << "/" << batch.size()
    //             << "  avg_pt=" << (ok_cnt ? sum_pt/ok_cnt : 0.0)
    //             << "  time=" << sw.RealTime() << " s\n";
    // }
}