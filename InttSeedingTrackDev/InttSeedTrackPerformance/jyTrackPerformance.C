#include "src/JYInttSeedTrackPerformance.h"
#include "src/JYInttSeedTrackPerformance.cxx"

#include "src/JYInttSeedTracking.h"
#include "src/JYInttSeedTracking.cxx"
#include "src/SPHTracKuma.h"

#include <typeinfo>

void jyTrackPerformance(Int_t runNum=0)
{
    // work area
    std::string fDir = "/mnt/e/sphenix/INTT-EMCAL/InttSeedingTrackDev/ParticleGen/output";

    // input file
    TChain *tc = new TChain("tree");  
    std::string fInputName = fDir + "/singleE1500to2500MeVEta0.root";
    tc->Add(fInputName.c_str()); 

    // output file
    std::string fOutputName = fDir+"/trackingWInttCalClu1500to2500MeVEtaWide_v2"; 

    TTree *tt = (TTree*)tc;
    InttSeedTrackPerformance *h = new InttSeedTrackPerformance(tt, fInputName, fOutputName, runNum);

    h->Loop(runNum);

    return;
}