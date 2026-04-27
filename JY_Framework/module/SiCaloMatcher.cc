#include "SiCaloMatcher.h"

#include "SiCaloMatchedTrack.h"
#include "SiCaloMatchedTrackContainer.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/getClass.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHIODataNode.h>
#include <phool/PHNodeIterator.h>
#include <phool/PHObject.h>
#include <phool/phool.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrCluster.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrDefs.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxTrackState.h>
#include <trackbase_historic/TrackSeed.h>

#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
#include <calobase/RawTowerGeomContainer.h>

#include <cmath>
#include <iostream>
#include <limits>

SiCaloMatcher::SiCaloMatcher(const std::string &name)
  : SubsysReco(name)
{
}

int SiCaloMatcher::Init(PHCompositeNode *topNode)
{
    return createOutputNode(topNode);
}

int SiCaloMatcher::InitRun(PHCompositeNode *topNode)
{
    return getNodes(topNode);
}

int SiCaloMatcher::End(PHCompositeNode * /*topNode*/)
{
    return Fun4AllReturnCodes::EVENT_OK;
}

int SiCaloMatcher::createOutputNode(PHCompositeNode *topNode)
{
    PHNodeIterator iter(topNode);

    PHCompositeNode *dstNode =
        dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));

    if (!dstNode)
    {
        std::cerr << PHWHERE << "DST node missing" << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    PHNodeIterator dstIter(dstNode);

    PHCompositeNode *svtxNode =
        dynamic_cast<PHCompositeNode *>(dstIter.findFirst("PHCompositeNode", "SVTX"));

    if (!svtxNode)
    {
        svtxNode = new PHCompositeNode("SVTX");
        dstNode->addNode(svtxNode);
    }

    m_outputContainer =
        findNode::getClass<SiCaloMatchedTrackContainer>(topNode, m_outputNodeName);

    if (!m_outputContainer)
    {
        m_outputContainer = new SiCaloMatchedTrackContainer();

        PHIODataNode<PHObject> *node =
            new PHIODataNode<PHObject>(
                m_outputContainer,
                m_outputNodeName,
                "PHObject");

        svtxNode->addNode(node);

        if (Verbosity() > 0)
        {
            std::cout << PHWHERE
                      << " created node SVTX/" << m_outputNodeName
                      << std::endl;
        }
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int SiCaloMatcher::getNodes(PHCompositeNode *topNode)
{
    m_trackMap = findNode::getClass<SvtxTrackMap>(topNode, m_trackMapName);
    if (!m_trackMap)
    {
        std::cout << PHWHERE << " cannot find " << m_trackMapName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    m_emcalClusters =
        findNode::getClass<RawClusterContainer>(topNode, m_emcalClusterName);

    if (!m_emcalClusters)
    {
        std::cout << PHWHERE << " cannot find " << m_emcalClusterName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    if (!m_use_emcal_radius)
    {
        m_emcalGeom =
            findNode::getClass<RawTowerGeomContainer>(topNode, m_emcalGeomName);

        if (!m_emcalGeom)
        {
            std::cout << PHWHERE << " cannot find " << m_emcalGeomName << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
    }

    if (m_store_silicon_hits)
    {
        m_trkrClusters =
            findNode::getClass<TrkrClusterContainer>(topNode, m_clusterMapName);

        if (!m_trkrClusters)
        {
            std::cout << PHWHERE << " cannot find " << m_clusterMapName << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }

        m_actsGeometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

        if (!m_actsGeometry)
        {
            std::cout << PHWHERE << " cannot find ActsGeometry" << std::endl;
            return Fun4AllReturnCodes::ABORTEVENT;
        }
    }

    m_outputContainer =
        findNode::getClass<SiCaloMatchedTrackContainer>(topNode, m_outputNodeName);

    if (!m_outputContainer)
    {
        std::cout << PHWHERE << " cannot find " << m_outputNodeName << std::endl;
        return Fun4AllReturnCodes::ABORTEVENT;
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

int SiCaloMatcher::process_event(PHCompositeNode *topNode)
{
    const int ret = getNodes(topNode);
    if (ret != Fun4AllReturnCodes::EVENT_OK)
    {
        return ret;
    }

    m_outputContainer->Reset();

    float caloRadiusEMCal = m_emcal_radius_user;
    if (!m_use_emcal_radius && m_emcalGeom)
    {
        caloRadiusEMCal = m_emcalGeom->get_radius();
    }

    unsigned int n_inserted = 0;

    for (auto &trackIter : *m_trackMap)
    {
        SvtxTrack *track = trackIter.second;

        if (!checkTrack(track))
        {
            continue;
        }

        SvtxTrackState *state = track->get_state(caloRadiusEMCal);

        if (!state)
        {
            if (Verbosity() > 1)
            {
                std::cout << PHWHERE
                          << " no projected state at R = "
                          << caloRadiusEMCal
                          << " for track " << track->get_id()
                          << std::endl;
            }
            continue;
        }

        const float proj_x = state->get_x();
        const float proj_y = state->get_y();
        const float proj_z = state->get_z();

        const float proj_phi = std::atan2(proj_y, proj_x);
        const float proj_eta = getEta(proj_x, proj_y, proj_z);

        RawCluster *best_cluster = nullptr;
        unsigned int best_cluster_id =
            std::numeric_limits<unsigned int>::max();

        float best_dphi = NAN;
        float best_dz = NAN;
        float best_deta = NAN;
        float best_dr = NAN;
        float best_score = std::numeric_limits<float>::max();

        RawClusterContainer::Range range = m_emcalClusters->getClusters();

        for (RawClusterContainer::Iterator clusIter = range.first;
             clusIter != range.second;
             ++clusIter)
        {
            RawCluster *cluster = clusIter->second;

            if (!cluster)
            {
                continue;
            }

            if (cluster->get_energy() < m_emcal_e_low_cut)
            {
                continue;
            }

            const float cx = cluster->get_x();
            const float cy = cluster->get_y();
            const float cz = cluster->get_z();

            const float cr = std::sqrt(cx * cx + cy * cy);
            if (cr <= 0)
            {
                continue;
            }

            const float cphi = std::atan2(cy, cx);
            const float ceta = getEta(cx, cy, cz);

            // Same convention as TrackCaloMatch:
            // scale cluster z to the chosen EMCal projection radius.
            const float radius_scale = caloRadiusEMCal / cr;
            const float cz_scaled = radius_scale * cz;

            const float dphi = piRange(proj_phi - cphi);
            const float dz = proj_z - cz_scaled;
            const float deta = proj_eta - ceta;
            const float dr = std::sqrt(dphi * dphi + deta * deta);

            if (std::fabs(dphi) > m_dphi_cut)
            {
                continue;
            }

            if (std::fabs(dz) > m_dz_cut)
            {
                continue;
            }

            const float score =
                std::sqrt((dphi / m_dphi_cut) * (dphi / m_dphi_cut) +
                          (dz / m_dz_cut) * (dz / m_dz_cut));

            if (score < best_score)
            {
                best_score = score;
                best_cluster = cluster;
                best_cluster_id = clusIter->first;
                best_dphi = dphi;
                best_dz = dz;
                best_deta = deta;
                best_dr = dr;
            }
        }

        if (!best_cluster)
        {
            continue;
        }

        SiCaloMatchedTrack out;
        out.Reset();

        out.track_id = track->get_id();
        out.charge = track->get_charge();

        out.track_pt = track->get_pt();
        out.track_eta = track->get_eta();
        out.track_phi = track->get_phi();
        out.track_px = track->get_px();
        out.track_py = track->get_py();
        out.track_pz = track->get_pz();

        out.track_x = track->get_x();
        out.track_y = track->get_y();
        out.track_z = track->get_z();

        out.track_quality = track->get_quality();
        out.track_crossing = track->get_crossing();

        out.proj_emcal_x = proj_x;
        out.proj_emcal_y = proj_y;
        out.proj_emcal_z = proj_z;
        out.proj_emcal_phi = proj_phi;
        out.proj_emcal_eta = proj_eta;

        out.has_emcal_match = true;
        out.emcal_cluster_id = best_cluster_id;

        out.emcal_x = best_cluster->get_x();
        out.emcal_y = best_cluster->get_y();
        out.emcal_z = best_cluster->get_z();
        out.emcal_phi = std::atan2(out.emcal_y, out.emcal_x);
        out.emcal_eta = getEta(out.emcal_x, out.emcal_y, out.emcal_z);

        out.emcal_e = best_cluster->get_energy();
        out.emcal_ecore = best_cluster->get_ecore();
        out.emcal_chi2 = best_cluster->get_chi2();
        out.emcal_prob = best_cluster->get_prob();

        out.match_dphi = best_dphi;
        out.match_dz = best_dz;
        out.match_deta = best_deta;
        out.match_dr = best_dr;
        out.match_score = best_score;

        if (m_store_silicon_hits)
        {
            fillSiliconHits(track, out);
        }

        m_outputContainer->insert(out.track_id, out);
        ++n_inserted;
    }

    if (Verbosity() > 0)
    {
        std::cout << PHWHERE
                  << " inserted " << n_inserted
                  << " matched SiCalo objects into "
                  << m_outputNodeName
                  << std::endl;
    }

    return Fun4AllReturnCodes::EVENT_OK;
}

bool SiCaloMatcher::checkTrack(SvtxTrack *track) const
{
    if (!track)
    {
        return false;
    }

    if (track->get_pt() < m_track_pt_low_cut)
    {
        return false;
    }

    return true;
}

float SiCaloMatcher::piRange(float dphi) const
{
    if (dphi > M_PI)
    {
        dphi -= 2.0 * M_PI;
    }

    if (dphi < -M_PI)
    {
        dphi += 2.0 * M_PI;
    }

    return dphi;
}

float SiCaloMatcher::getEta(float x, float y, float z) const
{
    const float r = std::sqrt(x * x + y * y);

    if (r <= 0)
    {
        return std::numeric_limits<float>::quiet_NaN();
    }

    return std::asinh(z / r);
}

void SiCaloMatcher::fillSiliconHits(
    SvtxTrack *track,
    SiCaloMatchedTrack &out)
{
    if (!track || !m_trkrClusters || !m_actsGeometry)
    {
        return;
    }

    TrackSeed *si_seed = track->get_silicon_seed();

    if (!si_seed)
    {
        if (Verbosity() > 1)
        {
            std::cout << PHWHERE
                      << " track " << track->get_id()
                      << " has no silicon seed"
                      << std::endl;
        }
        return;
    }

    for (auto keyIter = si_seed->begin_cluster_keys();
         keyIter != si_seed->end_cluster_keys();
         ++keyIter)
    {
        const TrkrDefs::cluskey key = *keyIter;

        TrkrCluster *cluster = m_trkrClusters->findCluster(key);

        if (!cluster)
        {
            continue;
        }

        Acts::Vector3 global =
            m_actsGeometry->getGlobalPosition(key, cluster);

        out.si_hit_key.push_back(static_cast<unsigned long long>(key));
        out.si_hit_detid.push_back(static_cast<int>(TrkrDefs::getTrkrId(key)));
        out.si_hit_layer.push_back(static_cast<int>(TrkrDefs::getLayer(key)));
        out.si_hit_x.push_back(global.x());
        out.si_hit_y.push_back(global.y());
        out.si_hit_z.push_back(global.z());
    }
}