#include "sbncode/SinglePhotonAnalysis/HelperFunctions/helper_gadget.h"
//table printers
std::vector<int> single_photon::Printer_header( std::vector< std::string> headings){

  std::vector<int> spacers;
  for( size_t index = 0; index < headings.size(); index++){
    std::cout<<headings[index];
    spacers.push_back( headings[index].size());
  }
  std::cout<<std::endl;
  return spacers;
}

void single_photon::Printer_content( std::vector< std::string > nums, std::vector<int> spacers){

  if(nums.size() != spacers.size()) {
    std::cout<<"CANNOT PRINT!"<<std::endl;
    return;
  }

  for( size_t index = 0; index < nums.size(); index++){
    std::cout<<std::setw(spacers[index])<<nums[index];
  }
  std::cout<<std::endl;
}
