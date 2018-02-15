#script to generate skeleton function for calculating weight
#usage:
# generate_weightcalc.py FUNCNAME
#
#output: 
# script creates FUNCNAMEWeightCalc.cxx
#
#to run your code, put it in generated cxx file
#and add the config to fcl file:

import os
import sys, getopt

def main(argv):
   funcname = ''

   try:
      opts, args = getopt.getopt(argv,"hf:",["funcname="])
   except getopt.GetoptError:
      print 'generate_weightcalc.py -f FUNCNAME'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print 'generate_weightcalc.py -f FUNCNAME'
         sys.exit()
      elif opt == '-f':
         funcname = arg
      else:
         print 'option not recognized'
         sys.exit(2)

   if ( funcname=='' ):
       print "You must specify function name with option -f!"
       sys.exit(2)

   filename="".join(x for x in funcname if x.isalnum())
   filename=filename+"WeightCalc.cxx"

   if (os.path.isfile(filename)):
       print "File",filename,"already exists!"
       sys.exit(2)

   print 'Genarating function', filename

   ofstr='''
#include \"WeightCalcCreator.h\"
#include \"WeightCalc.h\"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

#include "CLHEP/Random/RandGaussQ.h"

namespace sbncode {
  namespace evwgh {

class %(funcname)sWeightCalc : public WeightCalc {
   public:
     %(funcname)sWeightCalc();
     void Configure(fhicl::ParameterSet const& p);
     std::vector<std::vector<double> > GetWeight(art::Event & e);
   private:
     CLHEP::RandGaussQ *fGaussRandom;
     
   DECLARE_WEIGHTCALC(%(funcname)sWeightCalc)
};
%(funcname)sWeightCalc::%(funcname)sWeightCalc() {}

void %(funcname)sWeightCalc::Configure(fhicl::ParameterSet const& p) {
  // Get configuration for this function
  fhicl::ParameterSet const &pset=p.get<fhicl::ParameterSet> (GetName());

  // Prepare random generator
  art::ServiceHandle<art::RandomNumberGenerator> rng;
  fGaussRandom = new CLHEP::RandGaussQ(rng->getEngine(GetName()));    
}

std::vector<std::vector<double> > %(funcname)sWeightCalc::GetWeight(art::Event & e) {
  // Calculate weight(s) here
  std::vector<std::vector<double> > weight;
  return weight;
}

REGISTER_WEIGHTCALC(%(funcname)sWeightCalc)

  }
}
'''%{'funcname':funcname}

   of=open(filename,'w')

   of.write(ofstr)
   of.close()

if __name__ == "__main__":
   main(sys.argv[1:])

