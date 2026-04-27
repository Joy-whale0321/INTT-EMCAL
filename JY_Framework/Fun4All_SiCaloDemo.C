// Fun4All_SiCaloMatch.C

#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>

#include <G4_ActsGeom.C>

#include <trackreco/PHActsTrackProjection.h>
#include <trackbase_historic/SvtxTrack.h>

#include <SiCaloMatcher.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libg4dst.so)
R__LOAD_LIBRARY(libtrack_io.so)
R__LOAD_LIBRARY(libtrack_reco.so)
R__LOAD_LIBRARY(libcalo_io.so)
R__LOAD_LIBRARY(libcalo_reco.so)

// Change this to your compiled library name.
R__LOAD_LIBRARY(libSiCaloMatch.so)

int Fun4All_SiCaloMatch(
    const int nEvents = 100,
    const std::string &input_track_dst = "track_dst.list",
    const std::string &input_calo_dst = "calo_dst.list",
    const std::string &output_dst = "SiCaloMatch_nanoDST.root")
{
    Fun4AllServer *se = Fun4AllServer::instance();
    se->Verbosity(1);

    // -------------------------------------------------------
    // Geometry
    // -------------------------------------------------------
    G4Init();
    G4_ActsGeomInit();

    // -------------------------------------------------------
    // Input: track DST
    // -------------------------------------------------------
    Fun4AllInputManager *track_in = new Fun4AllDstInputManager("TRACK_IN");
    track_in->AddListFile(input_track_dst);
    se->registerInputManager(track_in);

    // -------------------------------------------------------
    // Input: calo DST
    // -------------------------------------------------------
    Fun4AllInputManager *calo_in = new Fun4AllDstInputManager("CALO_IN");
    calo_in->AddListFile(input_calo_dst);
    se->registerInputManager(calo_in);

    // -------------------------------------------------------
    // Projection to EMCal radius
    //
    // SiCaloMatcher will read:
    //     track->get_state(emcal_radius)
    // -------------------------------------------------------
    const float emcal_radius = 100.0;  // cm, tune by EMCal R scan

    auto projection = new PHActsTrackProjection("CaloProjection");
    projection->setLayerRadius(SvtxTrack::CEMC, emcal_radius);
    se->registerSubsystem(projection);

    // -------------------------------------------------------
    // SiCalo match module
    // -------------------------------------------------------
    SiCaloMatcher *matcher = new SiCaloMatcher("SiCaloMatcher");

    matcher->setTrackMapName("SvtxTrackMap");

    // Choose the one that exists in your DST:
    //     CLUSTERINFO_CEMC
    //     TOPOCLUSTER_EMCAL
    matcher->setEMCalClusterName("CLUSTERINFO_CEMC");
    // matcher->setEMCalClusterName("TOPOCLUSTER_EMCAL");

    matcher->setClusterMapName("TRKR_CLUSTER");
    matcher->setOutputNodeName("SiCaloMatchedTrackContainer");

    matcher->useEMCalRadius(true);
    matcher->setEMCalRadius(emcal_radius);

    matcher->setTrackPtLowCut(0.2);
    matcher->setEMCalELowCut(0.1);
    matcher->setDPhiCut(0.15);
    matcher->setDZCut(10.0);

    matcher->setStoreSiliconHits(true);

    se->registerSubsystem(matcher);

    // -------------------------------------------------------
    // Output nano DST
    // -------------------------------------------------------
    Fun4AllDstOutputManager *out =
        new Fun4AllDstOutputManager("SICALO_MATCH_DST", output_dst);

    out->AddNode("SiCaloMatchedTrackContainer");

    // For debugging, you may also keep original nodes.
    // Comment these out later if you want a smaller nano DST.
    out->AddNode("SvtxTrackMap");
    out->AddNode("TRKR_CLUSTER");

    se->registerOutputManager(out);

    // -------------------------------------------------------
    // Run
    // -------------------------------------------------------
    se->run(nEvents);
    se->End();

    delete se;
    gSystem->Exit(0);

    return 0;
}