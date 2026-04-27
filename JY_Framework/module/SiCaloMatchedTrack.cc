#include "SiCaloMatchedTrack.h"

void SiCaloMatchedTrack::Reset()
{
    track_id = std::numeric_limits<unsigned int>::max();
    charge = 0;

    track_pt = NAN;
    track_eta = NAN;
    track_phi = NAN;
    track_px = NAN;
    track_py = NAN;
    track_pz = NAN;

    track_x = NAN;
    track_y = NAN;
    track_z = NAN;

    track_quality = NAN;
    track_crossing = 0;

    proj_emcal_x = NAN;
    proj_emcal_y = NAN;
    proj_emcal_z = NAN;
    proj_emcal_phi = NAN;
    proj_emcal_eta = NAN;

    has_emcal_match = false;
    emcal_cluster_id = std::numeric_limits<unsigned int>::max();

    emcal_x = NAN;
    emcal_y = NAN;
    emcal_z = NAN;
    emcal_phi = NAN;
    emcal_eta = NAN;

    emcal_e = NAN;
    emcal_ecore = NAN;
    emcal_chi2 = NAN;
    emcal_prob = NAN;

    match_dphi = NAN;
    match_dz = NAN;
    match_deta = NAN;
    match_dr = NAN;
    match_score = NAN;

    si_hit_key.clear();
    si_hit_detid.clear();
    si_hit_layer.clear();
    si_hit_x.clear();
    si_hit_y.clear();
    si_hit_z.clear();
}

int SiCaloMatchedTrack::isValid() const
{
    return track_id != std::numeric_limits<unsigned int>::max();
}

void SiCaloMatchedTrack::identify(std::ostream &os) const
{
    os << "SiCaloMatchedTrack:"
       << " track_id=" << track_id
       << " has_emcal_match=" << has_emcal_match
       << " emcal_cluster_id=" << emcal_cluster_id
       << " track_pt=" << track_pt
       << " emcal_e=" << emcal_e
       << " dphi=" << match_dphi
       << " dz=" << match_dz
       << " n_si_hits=" << si_hit_x.size()
       << std::endl;
}