Input DST
  ├─ Silicon-track SvtxTrackMap
  ├─ EMCal/HCal cluster container
  └─ TRKR_CLUSTER / ActsGeometry

        ↓

PHActsTrackProjection
  - 给 SvtxTrack 插入 EMCal/HCal 半径处的 track state

        ↓

SiCaloMatcher
  - 读 SvtxTrackMap
  - 读 calo cluster container
  - 做 track–cluster matching
  - 输出 SiCaloMatchContainer

        ↓

SiCaloMomentumReco
  - 读 SiCaloMatchContainer
  - 计算 pT_reco
  - 更新 px/py/pz/E
  - 输出 SiCaloRecoContainer 或更新同一个 container

        ↓

DST / ROOT Tree output