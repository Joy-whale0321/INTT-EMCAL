# SiCalo Match Module/Class Design Summary

This note summarizes the current **SiCalo matching-only design**.  
The goal of this stage is to read silicon-track information and calorimeter-cluster information from DST, perform track--EMCal matching, and save the matched result into a lightweight nano-DST-style Fun4All node for later momentum reconstruction.

Current workflow:

```text
Track DST + Calo DST
        ↓
PHActsTrackProjection
        ↓
SiCaloMatcher
        ↓
SiCaloMatchedTrackContainer
        ↓
SiCaloMatch_nanoDST.root
```

At this stage, **momentum update / pT reconstruction is not done yet**.  
That will be handled later by a separate module, for example `SiCaloPtReco`, which will read `SiCaloMatchedTrackContainer` and call `PtCalculator`.

---

## 1. Overview of the three core classes

The matching part uses three core classes:

```text
SiCaloMatchedTrack
SiCaloMatchedTrackContainer
SiCaloMatcher
```

Their relationship is:

```text
SiCaloMatcher
    reads:
        SvtxTrackMap
        RawClusterContainer
        TRKR_CLUSTER
        ActsGeometry
    creates:# SiCalo Match Module/Class Design Summary

This note summarizes the current **SiCalo matching-only design**.  
The goal of this stage is to read silicon-track information and calorimeter-cluster information from DST, perform track--EMCal matching, and save the matched result into a lightweight nano-DST-style Fun4All node for later momentum reconstruction.

Current workflow:

```text
Track DST + Calo DST
        ↓
PHActsTrackProjection
        ↓
SiCaloMatcher
        ↓
SiCaloMatchedTrackContainer
        ↓
SiCaloMatch_nanoDST.root
```

At this stage, **momentum update / pT reconstruction is not done yet**.  
That will be handled later by a separate module, for example `SiCaloPtReco`, which will read `SiCaloMatchedTrackContainer` and call `PtCalculator`.

---

## 1. Overview of the three core classes

The matching part uses three core classes:

```text
SiCaloMatchedTrack
SiCaloMatchedTrackContainer
SiCaloMatcher
```

Their relationship is:

```text
SiCaloMatcher
    reads:
        SvtxTrackMap
        RawClusterContainer
        TRKR_CLUSTER
        ActsGeometry
    creates:
        many SiCaloMatchedTrack objects
    stores them in:
        SiCaloMatchedTrackContainer
```

A useful analogy with standard sPHENIX objects is:

```text
SvtxTrack                  = one reconstructed track
SvtxTrackMap               = one event's collection of tracks

SiCaloMatchedTrack         = one matched silicon-track + calo-cluster object
SiCaloMatchedTrackContainer = one event's collection of SiCalo matched objects

SiCaloMatcher              = the Fun4All module that builds the matched objects
```

---

## 2. `SiCaloMatchedTrack`

### Purpose

`SiCaloMatchedTrack` is a **data object**. It represents one matched Si+Calo object:

```text
one silicon track
    +
one matched EMCal cluster
    +
track projection information
    +
matching residuals
    +
silicon hit positions
```

It is not an algorithm module. It only stores information needed by later modules.

---

### Main stored information

#### 2.1 Original track information

These variables are copied from `SvtxTrack`:

```cpp
unsigned int track_id;
int charge;

float track_pt;
float track_eta;
float track_phi;
float track_px;
float track_py;
float track_pz;

float track_x;
float track_y;
float track_z;

float track_quality;
short track_crossing;
```

These are needed because later pT reconstruction must know the original track direction and charge.

---

#### 2.2 Track projection at EMCal radius

These variables store the track state at the EMCal matching surface:

```cpp
float proj_emcal_x;
float proj_emcal_y;
float proj_emcal_z;
float proj_emcal_phi;
float proj_emcal_eta;
```

They come from:

```cpp
track->get_state(emcal_radius)
```

This state is expected to be produced by `PHActsTrackProjection` before `SiCaloMatcher` runs.

---

#### 2.3 Matched EMCal cluster information

These variables describe the selected best-matched EMCal cluster:

```cpp
bool has_emcal_match;
unsigned int emcal_cluster_id;

float emcal_x;
float emcal_y;
float emcal_z;
float emcal_phi;
float emcal_eta;

float emcal_e;
float emcal_ecore;
float emcal_chi2;
float emcal_prob;
```

The variables are copied from the corresponding `RawCluster` in the EMCal cluster container.

---

#### 2.4 Matching residuals

These variables quantify the spatial difference between the projected track and the matched cluster:

```cpp
float match_dphi;
float match_dz;
float match_deta;
float match_dr;
float match_score;
```

Current matching uses mainly:

```text
|dphi| < dphi_cut
|dz|   < dz_cut
```

If multiple clusters pass the cuts, the module chooses the one with the smallest normalized score:

```text
score = sqrt((dphi / dphi_cut)^2 + (dz / dz_cut)^2)
```

---

#### 2.5 Silicon hit information

These vectors store silicon hit positions associated with the track:

```cpp
std::vector<unsigned long long> si_hit_key;
std::vector<int> si_hit_detid;
std::vector<int> si_hit_layer;
std::vector<float> si_hit_x;
std::vector<float> si_hit_y;
std::vector<float> si_hit_z;
```

These are filled by following the track's silicon seed:

```text
SvtxTrack
    ↓
track->get_silicon_seed()
    ↓
cluster keys
    ↓
TRKR_CLUSTER
    ↓
ActsGeometry::getGlobalPosition(key, cluster)
    ↓
global silicon hit x/y/z
```

The reason for storing silicon hit positions directly is that later modules, such as pT reconstruction or ML-based reconstruction, can use the matched object without needing to re-read `TRKR_CLUSTER` and `ActsGeometry`.

---

### Internal methods

#### `Reset()`

Resets all scalar values to invalid defaults such as `NAN` or `max unsigned int`, and clears all silicon-hit vectors.

#### `isValid()`

Checks whether the object has a valid `track_id`.

#### `identify()`

Prints a short summary of the object for debugging.

---

## 3. `SiCaloMatchedTrackContainer`

### Purpose

`SiCaloMatchedTrackContainer` stores all `SiCaloMatchedTrack` objects in one event.

If one event has five tracks matched to EMCal clusters, then the container stores five `SiCaloMatchedTrack` objects.

---

### Internal data structure

Internally, it uses a map:

```cpp
std::map<unsigned int, SiCaloMatchedTrack> m_map;
```

The key is currently the `track_id`:

```text
m_map[track_id] = SiCaloMatchedTrack
```

This makes it easy for later modules to retrieve the matched object corresponding to a given track.

---

### Main operations

#### `insert(unsigned int key, const SiCaloMatchedTrack &obj)`

Inserts one matched object into the container. If the key already exists, the old object is replaced.

#### `get(unsigned int key)`

Returns the matched object associated with a given key, usually `track_id`.

#### `begin()` and `end()`

Allow loop access:

```cpp
for (auto iter = container->begin(); iter != container->end(); ++iter)
{
    auto &matched = iter->second;
}
```

#### `Reset()`

Clears the event-level map at the beginning of each new event.

#### `identify()`

Prints the size of the container for debugging.

---

## 4. `SiCaloMatcher`

### Purpose

`SiCaloMatcher` is the actual **Fun4All `SubsysReco` module** that performs matching.

It reads track and calo nodes from the DST, performs track--EMCal matching, fills `SiCaloMatchedTrack` objects, and stores them into `SiCaloMatchedTrackContainer`.

---

### Input nodes

The module reads:

```text
SvtxTrackMap
```

This contains the reconstructed silicon-track or silicon-seed-converted track objects.

```text
CLUSTERINFO_CEMC or TOPOCLUSTER_EMCAL
```

This contains EMCal clusters. The actual node name depends on the input DST production.

```text
TRKR_CLUSTER
```

This contains tracker clusters. It is needed only if silicon hit positions are stored.

```text
ActsGeometry
```

This is needed to convert tracker cluster keys and local cluster information into global x/y/z positions.

Optionally:

```text
TOWERGEOM_CEMC_DETAILED
```

This is only needed if the EMCal radius is taken from geometry rather than set manually.

---

### Output node

The module creates and fills:

```text
SVTX/SiCaloMatchedTrackContainer
```

This node can be written out by `Fun4AllDstOutputManager` into a nano-DST-style ROOT file.

---

### Main algorithm in `process_event()`

For each event, `SiCaloMatcher` does the following:

#### Step 1: Read nodes

```cpp
getNodes(topNode);
```

This obtains the track map, EMCal cluster container, tracker cluster container, geometry, and output container.

---

#### Step 2: Clear output container

```cpp
m_outputContainer->Reset();
```

This removes results from the previous event.

---

#### Step 3: Determine EMCal matching radius

Currently the macro sets:

```cpp
matcher->useEMCalRadius(true);
matcher->setEMCalRadius(100.0);
```

So the module uses:

```text
R_EMCal = 100 cm
```

This radius must be consistent with `PHActsTrackProjection`.

---

#### Step 4: Loop over tracks

```cpp
for (auto &trackIter : *m_trackMap)
```

Each track is checked with `checkTrack()`.

Currently the basic cut is:

```text
track->get_pt() > m_track_pt_low_cut
```

---

#### Step 5: Get projected track state at EMCal radius

The module reads:

```cpp
SvtxTrackState *state = track->get_state(caloRadiusEMCal);
```

This state is expected to exist because the macro runs:

```cpp
PHActsTrackProjection
```

before `SiCaloMatcher`.

The projected position is then converted into:

```text
proj_emcal_x, proj_emcal_y, proj_emcal_z
proj_emcal_phi, proj_emcal_eta
```

---

#### Step 6: Loop over EMCal clusters

The module loops over all EMCal clusters:

```cpp
RawClusterContainer::Range range = m_emcalClusters->getClusters();
```

Each cluster must pass the energy cut:

```text
cluster->get_energy() > m_emcal_e_low_cut
```

---

#### Step 7: Compute residuals

For each track--cluster pair, it computes:

```text
dphi = phi_projected_track - phi_cluster
```

with proper wrapping to `[-pi, pi]`.

It also computes:

```text
deta = eta_projected_track - eta_cluster

dr = sqrt(dphi^2 + deta^2)
```

For `dz`, the cluster z is scaled to the chosen EMCal matching radius, following the same convention used in `TrackCaloMatch`:

```text
cluster_z_scaled = cluster_z * (R_EMCal / cluster_r)

dz = z_projected_track - cluster_z_scaled
```

---

#### Step 8: Apply matching cuts

Current matching requirements are:

```text
|dphi| < m_dphi_cut
|dz|   < m_dz_cut
```

Example settings in the macro:

```cpp
matcher->setDPhiCut(0.15);
matcher->setDZCut(10.0);
```

---

#### Step 9: Select best matched cluster

If multiple clusters pass the cuts, the best cluster is selected using:

```text
score = sqrt((dphi / dphi_cut)^2 + (dz / dz_cut)^2)
```

The cluster with the smallest score is selected.

---

#### Step 10: Fill `SiCaloMatchedTrack`

For each matched track, the module creates:

```cpp
SiCaloMatchedTrack out;
```

It fills:

```text
track information
projected track information
matched EMCal cluster information
matching residuals
silicon hit positions
```

---

#### Step 11: Insert into output container

The matched object is stored using:

```cpp
m_outputContainer->insert(out.track_id, out);
```

The key is the track id.

---

## 5. Macro-level workflow

The matching macro performs this sequence:

```text
1. Read track DST
2. Read calo DST
3. Initialize geometry
4. Run PHActsTrackProjection
5. Run SiCaloMatcher
6. Write SiCaloMatchedTrackContainer into nano DST
```

The important ordering is:

```text
PHActsTrackProjection must run before SiCaloMatcher
```

because `SiCaloMatcher` expects to find the projected state through:

```cpp
track->get_state(emcal_radius)
```

---

## 6. Current design choice

The current design intentionally separates matching from momentum reconstruction:

```text
SiCaloMatcher
    ↓
SiCaloMatchedTrackContainer
    ↓
SiCaloPtReco, later
```

This is useful because pT reconstruction will use `PtCalculator`, which has multiple methods. Keeping pT reconstruction separate prevents the matching module from becoming too complicated.

The current matching module only prepares the input object needed by later pT reconstruction.

---

## 7. Meaning of nano DST in this workflow

Here, **nano DST** means a ROOT file written by `Fun4AllDstOutputManager` that stores Fun4All node objects, such as:

```text
SVTX/SiCaloMatchedTrackContainer
```

This is different from a plain analysis `TTree`.

A plain analysis ROOT tree is convenient for plotting and ML, but it is not a Fun4All reconstruction object.  
The nano DST is intended to be readable by a later Fun4All macro/module, such as a future `SiCaloPtReco` module.

---

## 8. Next step

The next module should be something like:

```text
SiCaloPtReco
```

It should read:

```text
SiCaloMatchedTrackContainer
```

Then it should use `PtCalculator` to compute one or more reconstructed pT values, for example:

```text
pt_emd
pt_eproj
pt_mlemd
pt_mleproj
pt_selected
```

Then it can either:

1. write a new reco container, or  
2. extend the existing matched object/container design to include reco momentum fields.

The first version of `SiCaloMatcher` does not do this yet.

        many SiCaloMatchedTrack objects
    stores them in:
        SiCaloMatchedTrackContainer
```

A useful analogy with standard sPHENIX objects is:

```text
SiCaloMatchedTrack         = one matched silicon-track + calo-cluster object
SiCaloMatchedTrackContainer = one event's collection of SiCalo matched objects

SiCaloMatcher              = the Fun4All module that builds the matched objects
```

---

## 2. `SiCaloMatchedTrack`

### Purpose

`SiCaloMatchedTrack` is a **data object**. It represents one matched Si+Calo object:

```text
one silicon track
    +
one matched EMCal cluster
    +
track projection information
    +
matching residuals
    +
silicon hit positions
```

It is not an algorithm module. It only stores information needed by later modules.

---

### Main stored information

#### 2.1 Original track information

These variables are copied from `SvtxTrack`:

```cpp
unsigned int track_id;
int charge;

float track_pt;
float track_eta;
float track_phi;
float track_px;
float track_py;
float track_pz;

float track_x;
float track_y;
float track_z;

float track_quality;
short track_crossing;
```

These are needed because later pT reconstruction must know the original track direction and charge.

---

#### 2.2 Track projection at EMCal radius

These variables store the track state at the EMCal matching surface:

```cpp
float proj_emcal_x;
float proj_emcal_y;
float proj_emcal_z;
float proj_emcal_phi;
float proj_emcal_eta;
```

They come from:

```cpp
track->get_state(emcal_radius)
```

This state is expected to be produced by `PHActsTrackProjection` before `SiCaloMatcher` runs.

---

#### 2.3 Matched EMCal cluster information

These variables describe the selected best-matched EMCal cluster:

```cpp
bool has_emcal_match;
unsigned int emcal_cluster_id;

float emcal_x;
float emcal_y;
float emcal_z;
float emcal_phi;
float emcal_eta;

float emcal_e;
float emcal_ecore;
float emcal_chi2;
float emcal_prob;
```

The variables are copied from the corresponding `RawCluster` in the EMCal cluster container.

---

#### 2.4 Matching residuals

These variables quantify the spatial difference between the projected track and the matched cluster:

```cpp
float match_dphi;
float match_dz;
float match_deta;
float match_dr;
float match_score;
```

Current matching uses mainly:

```text
|dphi| < dphi_cut
|dz|   < dz_cut
```

If multiple clusters pass the cuts, the module chooses the one with the smallest normalized score:

```text
score = sqrt((dphi / dphi_cut)^2 + (dz / dz_cut)^2)
```

---

#### 2.5 Silicon hit information

These vectors store silicon hit positions associated with the track:

```cpp
std::vector<unsigned long long> si_hit_key;
std::vector<int> si_hit_detid;
std::vector<int> si_hit_layer;
std::vector<float> si_hit_x;
std::vector<float> si_hit_y;
std::vector<float> si_hit_z;
```

These are filled by following the track's silicon seed:

```text
SvtxTrack
    ↓
track->get_silicon_seed()
    ↓
cluster keys
    ↓
TRKR_CLUSTER
    ↓
ActsGeometry::getGlobalPosition(key, cluster)
    ↓
global silicon hit x/y/z
```

The reason for storing silicon hit positions directly is that later modules, such as pT reconstruction or ML-based reconstruction, can use the matched object without needing to re-read `TRKR_CLUSTER` and `ActsGeometry`.

---

### Internal methods

#### `Reset()`

Resets all scalar values to invalid defaults such as `NAN` or `max unsigned int`, and clears all silicon-hit vectors.

#### `isValid()`

Checks whether the object has a valid `track_id`.

#### `identify()`

Prints a short summary of the object for debugging.

---

## 3. `SiCaloMatchedTrackContainer`

### Purpose

`SiCaloMatchedTrackContainer` stores all `SiCaloMatchedTrack` objects in one event.

If one event has five tracks matched to EMCal clusters, then the container stores five `SiCaloMatchedTrack` objects.

---

### Internal data structure

Internally, it uses a map:

```cpp
std::map<unsigned int, SiCaloMatchedTrack> m_map;
```

The key is currently the `track_id`:

```text
m_map[track_id] = SiCaloMatchedTrack
```

This makes it easy for later modules to retrieve the matched object corresponding to a given track.

---

### Main operations

#### `insert(unsigned int key, const SiCaloMatchedTrack &obj)`

Inserts one matched object into the container. If the key already exists, the old object is replaced.

#### `get(unsigned int key)`

Returns the matched object associated with a given key, usually `track_id`.

#### `begin()` and `end()`

Allow loop access:

```cpp
for (auto iter = container->begin(); iter != container->end(); ++iter)
{
    auto &matched = iter->second;
}
```

#### `Reset()`

Clears the event-level map at the beginning of each new event.

#### `identify()`

Prints the size of the container for debugging.

---

## 4. `SiCaloMatcher`

### Purpose

`SiCaloMatcher` is the actual **Fun4All `SubsysReco` module** that performs matching.

It reads track and calo nodes from the DST, performs track--EMCal matching, fills `SiCaloMatchedTrack` objects, and stores them into `SiCaloMatchedTrackContainer`.

---

### Input nodes

The module reads:

```text
SvtxTrackMap
```

This contains the reconstructed silicon-track or silicon-seed-converted track objects.

```text
CLUSTERINFO_CEMC or TOPOCLUSTER_EMCAL
```

This contains EMCal clusters. The actual node name depends on the input DST production.

```text
TRKR_CLUSTER
```

This contains tracker clusters. It is needed only if silicon hit positions are stored.

```text
ActsGeometry
```

This is needed to convert tracker cluster keys and local cluster information into global x/y/z positions.

Optionally:

```text
TOWERGEOM_CEMC_DETAILED
```

This is only needed if the EMCal radius is taken from geometry rather than set manually.

---

### Output node

The module creates and fills:

```text
SVTX/SiCaloMatchedTrackContainer
```

This node can be written out by `Fun4AllDstOutputManager` into a nano-DST-style ROOT file.

---

### Main algorithm in `process_event()`

For each event, `SiCaloMatcher` does the following:

#### Step 1: Read nodes

```cpp
getNodes(topNode);
```

This obtains the track map, EMCal cluster container, tracker cluster container, geometry, and output container.

---

#### Step 2: Clear output container

```cpp
m_outputContainer->Reset();
```

This removes results from the previous event.

---

#### Step 3: Determine EMCal matching radius

Currently the macro sets:

```cpp
matcher->useEMCalRadius(true);
matcher->setEMCalRadius(100.0);
```

So the module uses:

```text
R_EMCal = 100 cm
```

This radius must be consistent with `PHActsTrackProjection`.

---

#### Step 4: Loop over tracks

```cpp
for (auto &trackIter : *m_trackMap)
```

Each track is checked with `checkTrack()`.

Currently the basic cut is:

```text
track->get_pt() > m_track_pt_low_cut
```

---

#### Step 5: Get projected track state at EMCal radius

The module reads:

```cpp
SvtxTrackState *state = track->get_state(caloRadiusEMCal);
```

This state is expected to exist because the macro runs:

```cpp
PHActsTrackProjection
```

before `SiCaloMatcher`.

The projected position is then converted into:

```text
proj_emcal_x, proj_emcal_y, proj_emcal_z
proj_emcal_phi, proj_emcal_eta
```

---

#### Step 6: Loop over EMCal clusters

The module loops over all EMCal clusters:

```cpp
RawClusterContainer::Range range = m_emcalClusters->getClusters();
```

Each cluster must pass the energy cut:

```text
cluster->get_energy() > m_emcal_e_low_cut
```

---

#### Step 7: Compute residuals

For each track--cluster pair, it computes:

```text
dphi = phi_projected_track - phi_cluster
```

with proper wrapping to `[-pi, pi]`.

It also computes:

```text
deta = eta_projected_track - eta_cluster

dr = sqrt(dphi^2 + deta^2)
```

For `dz`, the cluster z is scaled to the chosen EMCal matching radius, following the same convention used in `TrackCaloMatch`:

```text
cluster_z_scaled = cluster_z * (R_EMCal / cluster_r)

dz = z_projected_track - cluster_z_scaled
```

---

#### Step 8: Apply matching cuts

Current matching requirements are:

```text
|dphi| < m_dphi_cut
|dz|   < m_dz_cut
```

Example settings in the macro:

```cpp
matcher->setDPhiCut(0.15);
matcher->setDZCut(10.0);
```

---

#### Step 9: Select best matched cluster

If multiple clusters pass the cuts, the best cluster is selected using:

```text
score = sqrt((dphi / dphi_cut)^2 + (dz / dz_cut)^2)
```

The cluster with the smallest score is selected.

---

#### Step 10: Fill `SiCaloMatchedTrack`

For each matched track, the module creates:

```cpp
SiCaloMatchedTrack out;
```

It fills:

```text
track information
projected track information
matched EMCal cluster information
matching residuals
silicon hit positions
```

---

#### Step 11: Insert into output container

The matched object is stored using:

```cpp
m_outputContainer->insert(out.track_id, out);
```

The key is the track id.

---

## 5. Macro-level workflow

The matching macro performs this sequence:

```text
1. Read track DST
2. Read calo DST
3. Initialize geometry
4. Run PHActsTrackProjection
5. Run SiCaloMatcher
6. Write SiCaloMatchedTrackContainer into nano DST
```

The important ordering is:

```text
PHActsTrackProjection must run before SiCaloMatcher
```

because `SiCaloMatcher` expects to find the projected state through:

```cpp
track->get_state(emcal_radius)
```

---

## 6. Current design choice

The current design intentionally separates matching from momentum reconstruction:

```text
SiCaloMatcher
    ↓
SiCaloMatchedTrackContainer
    ↓
SiCaloPtReco, later
```

This is useful because pT reconstruction will use `PtCalculator`, which has multiple methods. Keeping pT reconstruction separate prevents the matching module from becoming too complicated.

The current matching module only prepares the input object needed by later pT reconstruction.

---

## 7. Meaning of nano DST in this workflow

Here, **nano DST** means a ROOT file written by `Fun4AllDstOutputManager` that stores Fun4All node objects, such as:

```text
SVTX/SiCaloMatchedTrackContainer
```

This is different from a plain analysis `TTree`.

A plain analysis ROOT tree is convenient for plotting and ML, but it is not a Fun4All reconstruction object.  
The nano DST is intended to be readable by a later Fun4All macro/module, such as a future `SiCaloPtReco` module.

---

## 8. Next step

The next module should be something like:

```text
SiCaloPtReco
```

It should read:

```text
SiCaloMatchedTrackContainer
```

Then it should use `PtCalculator` to compute one or more reconstructed pT values, for example:

```text
pt_emd
pt_eproj
pt_mlemd
pt_mleproj
pt_selected
```

Then it can either:

1. write a new reco container, or  
2. extend the existing matched object/container design to include reco momentum fields.

The first version of `SiCaloMatcher` does not do this yet.
