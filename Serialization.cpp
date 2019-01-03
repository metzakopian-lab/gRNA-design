#include "Serialization.h"
#include <iostream>


void modelSerialize(const GuideModel& gRNAs , const std::string& filename)
{
        std::cerr << "Dumping results to " << filename.c_str() << std::endl;
        std::ofstream out(filename.c_str());
  
        {
                boost::archive::binary_oarchive oa(out);
                oa << gRNAs;
        }
}



GuideModel modelDeserialize(const std::string& filename)
{
        GuideModel gRNAs;
        std::cerr << "Reading Model from " << filename << std::endl;
        
        std::ifstream in(filename.c_str());
        {
                boost::archive::binary_iarchive ia(in);
                ia >> gRNAs;
        }
        return gRNAs;
}
