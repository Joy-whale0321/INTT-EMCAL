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
    std::string fInputName = fDir + "/singleE1500to2500MeVEta0.root";
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
    auto h_x = new TH1F("h_x", "X position;X (cm);Counts", 100, -100, 100);
    auto h_y = new TH1F("h_y", "Y position;Y (cm);Counts", 100, -100, 100);
    auto h_z = new TH1F("h_z", "Z position;Z (cm);Counts", 100, -300, 300);
    auto h_R = new TH1F("h_R", "Radial distance;R (cm);Counts", 100, 0, 150);
    auto h_towerphi = new TH1F("h_towerphi", "Azimuthal angle;Phi (rad);Counts", 100, -TMath::Pi(), TMath::Pi());
    auto h_clustphi = new TH1F("h_clustphi", "Azimuthal angle;Phi (rad);Counts", 100, -TMath::Pi(), TMath::Pi());
    auto h_eta = new TH1F("h_eta", "Polar angle;Theta (rad);Counts", 100, 0, TMath::Pi());

    // 设置分支地址
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
    for (Long64_t i = 0; i < 10; i++)
    {
        caloinfotree->GetEntry(i);

        Double_t TotEMCalE = 0.;
        Double_t ModifEMCalPhi = 0.;

        // 遍历 tower vector 中的每个元素
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

            if((towerval_edep > Ecalo_threshold)&&(towerval_R<100))
            {
                TotEMCalE += towerval_edep;
                ModifEMCalPhi += towerval_edep*towerval_phi;

                h_towerphi->Fill(towerval_phi);
            }           
        }
        ModifEMCalPhi /= TotEMCalE;
        cout<<"i is: "<< i <<", towerval_phi is: "<< ModifEMCalPhi <<endl;

        // 遍历 cluster vector 中的每个元素
        for (size_t j = 0; j < caloClus_x->size(); j++)
        {          
            double clusterval_x = caloClus_x->at(j);
            double clusterval_y = caloClus_y->at(j);
            double clusterval_z = caloClus_z->at(j);
            double clusterval_R = caloClus_R->at(j);
            double clusterval_phi = caloClus_phi->at(j);
            double clusterval_edep = caloClus_edep->at(j);

            double Ecalo_threshold = 0.5;
            if((clusterval_edep > Ecalo_threshold)&&(clusterval_R<100))
            {
                h_clustphi->Fill(clusterval_phi);
                cout<<"i is: "<< i <<", clusterval_phi is: "<< clusterval_phi <<endl;
            }           
        }
    }

    // 保存直方图到文件
    TFile outfile(fOutputName.c_str(), "RECREATE");
    
    h_towerphi->Write();
    h_clustphi->Write();

    std::cout << "Histograms saved to " << fOutputName << std::endl;
}


// // 对emcalhit取加权平均(w-energy)，energy加上ihcal和ohcal
// void CalESumAndCorrPosi(tracKuma trk, std::vector<hitStruct> vEmcalHits,\
//     std::vector<hitStruct> vIHCalHits, std::vector<hitStruct> vOHCalHits)
// {
//     Double_t refCalPhi = trk.getHitPhi(6);
//     Double_t refCalTheta = trk.getHitTheta(6);
//     Double_t TotEMCalE = 0.;
//     Double_t ModifEMCalPhi = 0.;
//     Double_t ModifEMCalTheta = 0.;

//     // 加权平均算x,y,z
//     for(Int_t iEmcal = 0; iEmcal < vEmcalHits.size(); iEmcal++)
//     {
//         Double_t hitTheta = 2*atan(std::exp(-vEmcalHits.at(iEmcal).eta));
//         if((hitTheta < refCalTheta - TMath::Pi()/20)&&(refCalTheta + TMath::Pi()/20 < hitTheta)) continue;
//         if((vEmcalHits.at(iEmcal).phi < refCalPhi - TMath::Pi()/20)&&(refCalPhi + TMath::Pi()/20 < vEmcalHits.at(iEmcal).phi)) continue;
        
//         TotEMCalE += vEmcalHits.at(iEmcal).energy;
//         ModifEMCalPhi += vEmcalHits.at(iEmcal).energy*vEmcalHits.at(iEmcal).phi;
//         ModifEMCalTheta += vEmcalHits.at(iEmcal).energy*hitTheta;
//     }
//     ModifEMCalPhi /= TotEMCalE;
//     ModifEMCalTheta /= TotEMCalE;

//     // 加权平均算R,phi,theta
//     for(Int_t iEmcal = 0; iEmcal < vEmcalHits.size(); iEmcal++)
//     {
//         Double_t hitTheta = 2*atan(std::exp(-vEmcalHits.at(iEmcal).eta));
//         if((hitTheta < refCalTheta - TMath::Pi()/20)&&(refCalTheta + TMath::Pi()/20 < hitTheta)) continue;
//         if((vEmcalHits.at(iEmcal).phi < refCalPhi - TMath::Pi()/20)&&(refCalPhi + TMath::Pi()/20 < vEmcalHits.at(iEmcal).phi)) continue;
        
//         TotEMCalE += vEmcalHits.at(iEmcal).energy;
//         ModifEMCalPhi += vEmcalHits.at(iEmcal).energy*vEmcalHits.at(iEmcal).phi;
//         ModifEMCalTheta += vEmcalHits.at(iEmcal).energy*hitTheta;
//     }
//     ModifEMCalPhi /= TotEMCalE;
//     ModifEMCalTheta /= TotEMCalE;
// }