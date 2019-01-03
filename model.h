#pragma once 

#include <boost/serialization/base_object.hpp>
#include <map>
#include <vector>


class GuideMeta
{
  
  friend class boost::serialization::access;
 public:
  size_t start_pos;
  size_t end_pos;
  size_t id;
  std::string chromosome;

};


typedef std::map<std::string, std::vector<GuideMeta>> GuideModel;
