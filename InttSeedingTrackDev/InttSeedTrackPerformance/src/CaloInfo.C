#include <iostream>           
#include <string>             
#include <vector>
#include <TChain.h>           
#include <TTree.h>            
#include <TH1F.h>             
#include <TFile.h>            
#include <TMath.h>            

void CaloInfo()
{
    // 工作目录
    std::string fDir = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ParticleGen/output";

    // 输入文件
    TChain *tc = new TChain("tree");  
    std::string fInputName = fDir + "/ana457_e_1_10GeV.root";
    tc->Add(fInputName.c_str()); 

    // 输出文件
    std::string fOutputName = fDir + "/calo_position_histograms.root"; 

    TTree *caloinfotree = (TTree*)tc;
    if (!caloinfotree)
    {
        std::cerr << "Error: Input tree is null!" << std::endl;
        return;
    }

    // 定义变量
    std::vector<double> *CEMC_Hit_x = nullptr; 
    std::vector<double> *CEMC_Hit_y = nullptr; 
    std::vector<double> *CEMC_Hit_z = nullptr;
    std::vector<double> *CEMC_Hit_Edep = nullptr;

    std::vector<double> *tower_x = nullptr;
    std::vector<double> *tower_y = nullptr;
    std::vector<double> *tower_z = nullptr;
    std::vector<double> *tower_phi = nullptr;
    std::vector<double> *tower_eta = nullptr;
    std::vector<double> *tower_edep = nullptr;

    std::vector<double> *caloClus_x  = nullptr; 
    std::vector<double> *caloClus_y  = nullptr;
    std::vector<double> *caloClus_z  = nullptr;
    std::vector<double> *caloClus_R  = nullptr;
    std::vector<double> *caloClus_phi  = nullptr;
    std::vector<double> *caloClus_edep = nullptr;

    // 定义直方图
    auto h_g4hitphi = new TH1F("h_g4hitphi", "Azimuthal angle;Phi (rad);Counts", 100, -TMath::Pi(), TMath::Pi());
    auto h_towerphi = new TH1F("h_towerphi", "Azimuthal angle;Phi (rad);Counts", 100, -TMath::Pi(), TMath::Pi());
    auto h_clustphi = new TH1F("h_clustphi", "Azimuthal angle;Phi (rad);Counts", 100, -TMath::Pi(), TMath::Pi());
    auto h_deltatowertruth = new TH1F("h_deltatowertruth", "h_deltatowertruth;R;Counts", 100, 0, 20);
    auto h_deltaclusttruth = new TH1F("h_deltaclusttruth", "h_deltaclusttruth;R;Counts", 100, 0, 20);
    auto h_deltaclusttower = new TH1F("h_deltaclusttower", "h_deltaclusttower;R;Counts", 500, 0, 10);

    auto h2_zr_g4hit = new TH2D("h2_zr_g4hit", "g4hit energy-weighted RZ distribution;Z (cm);R (cm);Weighted Counts", 100, -150, 150, 200, 0, 200);
    auto h2_zr_tower = new TH2D("h2_zr_tower", "tower energy-weighted RZ distribution;Z (cm);R (cm);Weighted Counts", 100, -150, 150, 200, 0, 200);
    auto h2_zr_clust = new TH2D("h2_zr_clust", "clust energy-weighted RZ distribution;Z (cm);R (cm);Weighted Counts", 100, -150, 150, 200, 0, 200);



    // 设置分支地址
    caloinfotree->SetBranchAddress("CEMC_Hit_x", &CEMC_Hit_x);
    caloinfotree->SetBranchAddress("CEMC_Hit_y", &CEMC_Hit_y);
    caloinfotree->SetBranchAddress("CEMC_Hit_z", &CEMC_Hit_z);
    caloinfotree->SetBranchAddress("CEMC_Hit_Edep", &CEMC_Hit_Edep);

    caloinfotree->SetBranchAddress("tower_X", &tower_x);
    caloinfotree->SetBranchAddress("tower_Y", &tower_y);
    caloinfotree->SetBranchAddress("tower_Z", &tower_z);
    caloinfotree->SetBranchAddress("tower_Phi", &tower_phi);
    caloinfotree->SetBranchAddress("tower_Eta", &tower_eta);
    caloinfotree->SetBranchAddress("tower_edep", &tower_edep);

    caloinfotree->SetBranchAddress("caloClus_X", &caloClus_x);
    caloinfotree->SetBranchAddress("caloClus_Y", &caloClus_y);
    caloinfotree->SetBranchAddress("caloClus_Z", &caloClus_z);
    caloinfotree->SetBranchAddress("caloClus_R", &caloClus_R);
    caloinfotree->SetBranchAddress("caloClus_Phi", &caloClus_phi);
    caloinfotree->SetBranchAddress("caloClus_edep", &caloClus_edep);

    // 遍历树的条目
    Long64_t nentries = caloinfotree->GetEntries();
    // for (Long64_t i = 0; i < nentries; i++)
    for (Long64_t i = 0; i < nentries; i++)
    {
        caloinfotree->GetEntry(i);

        // 遍历 CEMC_Hit vector 中的每个元素
        Double_t g4hit_TotEMCalE = 0.;
        Double_t g4hit_ModifEMCalPhi = 0.;
        Double_t g4hit_ModifEMCal_x = 0.;
        Double_t g4hit_ModifEMCal_y = 0.;

        for (size_t j = 0; j < CEMC_Hit_x->size(); j++)
        {
            double hitval_x = CEMC_Hit_x->at(j);
            double hitval_y = CEMC_Hit_y->at(j);
            double hitval_z = CEMC_Hit_z->at(j);
            double hitval_R = TMath::Sqrt(hitval_x * hitval_x + hitval_y * hitval_y);
            double hitval_phi = TMath::ATan2(hitval_y, hitval_x);
            double hitval_eta = TMath::ATan2(hitval_R, hitval_z);
            double hitval_edep = CEMC_Hit_Edep->at(j);

            if(hitval_R<300)
            {
                g4hit_TotEMCalE += hitval_edep;
                g4hit_ModifEMCalPhi += hitval_edep * hitval_phi;

                g4hit_ModifEMCal_x += hitval_edep * hitval_x;
                g4hit_ModifEMCal_y += hitval_edep * hitval_y;

                h2_zr_g4hit->Fill(hitval_z, hitval_R, hitval_edep);
            }
        }
        g4hit_ModifEMCalPhi /= g4hit_TotEMCalE;
        g4hit_ModifEMCal_x /= g4hit_TotEMCalE;
        g4hit_ModifEMCal_y /= g4hit_TotEMCalE;
        double g4hit_ModifEMCal_phi = TMath::ATan2(g4hit_ModifEMCal_y, g4hit_ModifEMCal_x);

        h_g4hitphi->Fill(g4hit_ModifEMCalPhi);


        // 遍历 tower vector 中的每个元素
        Double_t tower_TotEMCalE = 0.;
        Double_t tower_ModifEMCalPhi = 0.;
        Double_t tower_ModifEMCal_x = 0.;
        Double_t tower_ModifEMCal_y = 0.;

        for (size_t j = 0; j < tower_x->size(); j++)
        {      
            double towerval_x = tower_x->at(j);
            double towerval_y = tower_y->at(j);
            double towerval_z = tower_z->at(j);
            double towerval_phi = tower_phi->at(j);
            double towerval_eta = tower_eta->at(j);
            double towerval_edep = tower_edep->at(j);
            
            double Ecalo_threshold = 0.1;
            double towerval_R = TMath::Sqrt(towerval_x * towerval_x + towerval_y * towerval_y);

            if((towerval_edep > Ecalo_threshold)&&(towerval_R<200))
            {
                tower_TotEMCalE += towerval_edep;

                tower_ModifEMCalPhi += towerval_edep * towerval_phi;
                tower_ModifEMCal_x += towerval_edep * towerval_x;
                tower_ModifEMCal_y += towerval_edep * towerval_y;

                h2_zr_tower->Fill(towerval_z, towerval_R, towerval_edep);
            }                        
        }
        tower_ModifEMCalPhi /= tower_TotEMCalE;
        tower_ModifEMCal_x /= tower_TotEMCalE;
        tower_ModifEMCal_y /= tower_TotEMCalE;
        double tower_ModifEMCal_phi = TMath::ATan2(tower_ModifEMCal_y, tower_ModifEMCal_x);

        h_towerphi->Fill(tower_ModifEMCalPhi);


        // 遍历 cluster vector 中的每个元素
        Double_t cluster_TotEMCalE = 0.;
        Double_t cluster_ModifEMCalPhi = 0.;
        Double_t cluster_ModifEMCal_x = 0.;
        Double_t cluster_ModifEMCal_y = 0.;

        for (size_t j = 0; j < caloClus_x->size(); j++)
        {          
            double clusterval_x = caloClus_x->at(j);
            double clusterval_y = caloClus_y->at(j);
            double clusterval_z = caloClus_z->at(j);
            double clusterval_R = caloClus_R->at(j);
            double clusterval_phi = caloClus_phi->at(j);
            double clusterval_edep = caloClus_edep->at(j);

            double Ecalo_threshold = 0.1;
            if((clusterval_edep > Ecalo_threshold)&&(clusterval_R<200))
            {
                h_clustphi->Fill(clusterval_phi);
                cluster_ModifEMCal_x += clusterval_x;
                cluster_ModifEMCal_y += clusterval_y;
                cluster_TotEMCalE += 1;

                h2_zr_clust->Fill(clusterval_z, clusterval_R, clusterval_edep);
            }           
        }
        cluster_ModifEMCal_x /= cluster_TotEMCalE;
        cluster_ModifEMCal_y /= cluster_TotEMCalE;

        // calculate delta truth tower cluster
        double delta_towertruth = sqrt((tower_ModifEMCal_x - g4hit_ModifEMCal_x)*(tower_ModifEMCal_x - g4hit_ModifEMCal_x) + (tower_ModifEMCal_y - g4hit_ModifEMCal_y)*(tower_ModifEMCal_y - g4hit_ModifEMCal_y)) ;
        double delta_clustertruth = sqrt((cluster_ModifEMCal_x - g4hit_ModifEMCal_x)*(cluster_ModifEMCal_x - g4hit_ModifEMCal_x) + (cluster_ModifEMCal_y - g4hit_ModifEMCal_y)*(cluster_ModifEMCal_y - g4hit_ModifEMCal_y)) ;
        
        double delta_clustertower = sqrt((cluster_ModifEMCal_x - tower_ModifEMCal_x)*(cluster_ModifEMCal_x - tower_ModifEMCal_x) + (cluster_ModifEMCal_y - tower_ModifEMCal_y)*(cluster_ModifEMCal_y - tower_ModifEMCal_y)) ;


        h_deltatowertruth->Fill(delta_towertruth);
        h_deltaclusttruth->Fill(delta_clustertruth);
        h_deltaclusttower->Fill(delta_clustertower);

    }

    // 保存直方图到文件
    TFile outfile(fOutputName.c_str(), "RECREATE");
    
    h_g4hitphi->Write();
    h_towerphi->Write();
    h_clustphi->Write();
    h_deltatowertruth->Write();
    h_deltaclusttruth->Write();
    h_deltaclusttower->Write();

    h2_zr_g4hit->Write();
    h2_zr_tower->Write();
    h2_zr_clust->Write();

    std::cout << "Histograms saved to " << fOutputName << std::endl;
}