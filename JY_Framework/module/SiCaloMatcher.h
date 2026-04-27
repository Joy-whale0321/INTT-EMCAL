#ifndef SICALOMATCHER_H
#define SICALOMATCHER_H

#include <fun4all/SubsysReco.h>

#include <string>

class PHCompositeNode;

class SvtxTrack;
class SvtxTrackMap;
class RawClusterContainer;
class RawTowerGeomContainer;
class TrkrClusterContainer;
class ActsGeometry;
class SiCaloMatchedTrack;
class SiCaloMatchedTrackContainer;

class SiCaloMatcher : public SubsysReco
{
 public:
    SiCaloMatcher(const std::string &name = "SiCaloMatcher");
    ~SiCaloMatcher() override = default;

    int Init(PHCompositeNode *topNode) override;
    int InitRun(PHCompositeNode *topNode) override;
    int process_event(PHCompositeNode *topNode) override;
    int End(PHCompositeNode *topNode) override;

    void setTrackMapName(const std::string &name) { m_trackMapName = name; }
    void setEMCalClusterName(const std::string &name) { m_emcalClusterName = name; }
    void setClusterMapName(const std::string &name) { m_clusterMapName = name; }
    void setEMCalGeomName(const std::string &name) { m_emcalGeomName = name; }
    void setOutputNodeName(const std::string &name) { m_outputNodeName = name; }

    void useEMCalRadius(bool use) { m_use_emcal_radius = use; }
    void setEMCalRadius(float r) { m_emcal_radius_user = r; }

    void setTrackPtLowCut(float pt) { m_track_pt_low_cut = pt; }
    void setEMCalELowCut(float e) { m_emcal_e_low_cut = e; }
    void setDPhiCut(float cut) { m_dphi_cut = cut; }
    void setDZCut(float cut) { m_dz_cut = cut; }

    void setStoreSiliconHits(bool value) { m_store_silicon_hits = value; }

 private:
    int createOutputNode(PHCompositeNode *topNode);
    int getNodes(PHCompositeNode *topNode);

    bool checkTrack(SvtxTrack *track) const;

    float piRange(float dphi) const;
    float getEta(float x, float y, float z) const;

    void fillSiliconHits(SvtxTrack *track, SiCaloMatchedTrack &out);

    std::string m_trackMapName = "SvtxTrackMap";
    std::string m_emcalClusterName = "CLUSTERINFO_CEMC";
    std::string m_clusterMapName = "TRKR_CLUSTER";
    std::string m_emcalGeomName = "TOWERGEOM_CEMC_DETAILED";
    std::string m_outputNodeName = "SiCaloMatchedTrackContainer";

    SvtxTrackMap *m_trackMap = nullptr;
    RawClusterContainer *m_emcalClusters = nullptr;
    RawTowerGeomContainer *m_emcalGeom = nullptr;
    TrkrClusterContainer *m_trkrClusters = nullptr;
    ActsGeometry *m_actsGeometry = nullptr;
    SiCaloMatchedTrackContainer *m_outputContainer = nullptr;

    bool m_use_emcal_radius = true;
    float m_emcal_radius_user = 100.0;

    float m_track_pt_low_cut = 0.2;
    float m_emcal_e_low_cut = 0.1;
    float m_dphi_cut = 0.15;
    float m_dz_cut = 10.0;

    bool m_store_silicon_hits = true;
};

#endif