#ifndef SICALOMATCHEDTRACKCONTAINER_H
#define SICALOMATCHEDTRACKCONTAINER_H

#include "SiCaloMatchedTrack.h"

#include <phool/PHObject.h>

#include <iostream>
#include <map>

class SiCaloMatchedTrackContainer : public PHObject
{
 public:
    using Map = std::map<unsigned int, SiCaloMatchedTrack>;
    using Iterator = Map::iterator;
    using ConstIterator = Map::const_iterator;

    SiCaloMatchedTrackContainer() = default;
    ~SiCaloMatchedTrackContainer() override = default;

    void Reset() override;
    int isValid() const override;
    void identify(std::ostream &os = std::cout) const override;

    SiCaloMatchedTrack *insert(unsigned int key, const SiCaloMatchedTrack &obj);
    SiCaloMatchedTrack *get(unsigned int key);
    const SiCaloMatchedTrack *get(unsigned int key) const;

    Iterator begin() { return m_map.begin(); }
    Iterator end() { return m_map.end(); }

    ConstIterator begin() const { return m_map.begin(); }
    ConstIterator end() const { return m_map.end(); }

    std::size_t size() const { return m_map.size(); }

 private:
    Map m_map;

    ClassDefOverride(SiCaloMatchedTrackContainer, 1)
};

#endif