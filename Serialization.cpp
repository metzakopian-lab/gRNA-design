#include "Serialization.h"
#include <iostream>


void stats(const GuideModel& gRNAs)
{
  std::cerr << "Found " << gRNAs.size() << " unique gRNAs" << std::endl;
  size_t max = 0, sum = 0;
  for (auto const& g : gRNAs)
  {
    
    max  = max > g.second.size() ? max : g.second.size();
    sum += g.second.size();
  }

  std::cerr << "Max duplicates: " << max << ", Total unique loci: " << sum <<std::endl;
}



void modelSerialize(const GuideModel& gRNAs , const std::string& filename)
{
  std::cerr << "Dumping results to " << filename.c_str() << std::endl;
  stats(gRNAs);
  std::ofstream out(filename.c_str());
  
  {
    boost::archive::binary_oarchive oa(out);
    oa << gRNAs;
  }
}




void modelDeserialize(GuideModel& gRNAs, const std::string& filename)
{

  std::cerr << "Reading Model from " << filename << std::endl;
  
  
  std::ifstream in(filename.c_str());
  {
    boost::archive::binary_iarchive ia(in);
    ia >> gRNAs;
  }
  stats(gRNAs);

}
