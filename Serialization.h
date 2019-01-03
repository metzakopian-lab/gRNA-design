#pragma once

#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <fstream>


#include "Model.h"


template <class Archive>
void serialize(Archive& ar, GuideMeta& obj, unsigned int version)
{
        ar & obj.start_pos;
        ar & obj.end_pos;
        ar & obj.id;
        ar & obj.chromosome;
    
}


void modelSerialize(const GuideModel& gRNAs , const std::string& filename);


void modelDeserialize(GuideModel& obj, const std::string& filename);

