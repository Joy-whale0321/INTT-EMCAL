#ifndef SICALOMATCHEDTRACK_H
#define SICALOMATCHEDTRACK_H

#include <phool/PHObject.h>

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

class SiCaloMatchedTrack : public PHObject
{
 public:
    SiCaloMatchedTrack() = default;
    ~SiCaloMatchedTrack() override = default;

    void Reset() override;
    int isValid() const override;
    void identify(std::ostream &os = std::cout) const override;

    // -------------------------------------------------------
    // Track identity and original track information
    // -------------------------------------------------------
    unsigned int track_id = std::numeric_limits<unsigned int>::max();
    int charge = 0;

    float track_pt = NAN;
    float track_eta = NAN;
    float track_phi = NAN;
    float track_px = NAN;
    float track_py = NAN;
    float track_pz = NAN;

    float track_x = NAN;
    float track_y = NAN;
    float track_z = NAN;

    float track_quality = NAN;
    short track_crossing = 0;

    // -------------------------------------------------------
    // Track projection at EMCal radius
    // -------------------------------------------------------
    float proj_emcal_x = NAN;
    float proj_emcal_y = NAN;
    float proj_emcal_z = NAN;
    float proj_emcal_phi = NAN;
    float proj_emcal_eta = NAN;

    // -------------------------------------------------------
    // Matched EMCal cluster information
    // -------------------------------------------------------
    bool has_emcal_match = false;
    unsigned int emcal_cluster_id = std::numeric_limits<unsigned int>::max();

    float emcal_x = NAN;
    float emcal_y = NAN;
    float emcal_z = NAN;
    float emcal_phi = NAN;
    float emcal_eta = NAN;

    float emcal_e = NAN;
    float emcal_ecore = NAN;
    float emcal_chi2 = NAN;
    float emcal_prob = NAN;

    // -------------------------------------------------------
    // Match residuals
    // -------------------------------------------------------
    float match_dphi = NAN;
    float match_dz = NAN;
    float match_deta = NAN;
    float match_dr = NAN;
    float match_score = NAN;

    // -------------------------------------------------------
    // Silicon hits associated with this track
    // Stored directly so next module does not need TRKR_CLUSTER.
    // -------------------------------------------------------
    std::vector<unsigned long long> si_hit_key;
    std::vector<int> si_hit_detid;
    std::vector<int> si_hit_layer;
    std::vector<float> si_hit_x;
    std::vector<float> si_hit_y;
    std::vector<float> si_hit_z;

    ClassDefOverride(SiCaloMatchedTrack, 1)
};

#endif