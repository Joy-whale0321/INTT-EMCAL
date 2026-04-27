#include "SiCaloMatchedTrackContainer.h"

void SiCaloMatchedTrackContainer::Reset()
{
    m_map.clear();
}

int SiCaloMatchedTrackContainer::isValid() const
{
    return 1;
}

void SiCaloMatchedTrackContainer::identify(std::ostream &os) const
{
    os << "SiCaloMatchedTrackContainer size = "
       << m_map.size()
       << std::endl;
}

SiCaloMatchedTrack *SiCaloMatchedTrackContainer::insert(
    unsigned int key,
    const SiCaloMatchedTrack &obj)
{
    auto ret = m_map.insert(std::make_pair(key, obj));

    if (!ret.second)
    {
        ret.first->second = obj;
    }

    return &(ret.first->second);
}

SiCaloMatchedTrack *SiCaloMatchedTrackContainer::get(unsigned int key)
{
    auto iter = m_map.find(key);

    if (iter == m_map.end())
    {
        return nullptr;
    }

    return &(iter->second);
}

const SiCaloMatchedTrack *SiCaloMatchedTrackContainer::get(unsigned int key) const
{
    auto iter = m_map.find(key);

    if (iter == m_map.end())
    {
        return nullptr;
    }

    return &(iter->second);
}