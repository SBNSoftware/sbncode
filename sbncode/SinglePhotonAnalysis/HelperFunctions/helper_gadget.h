#ifndef HELPER_GADGET_H
#define HELPER_GADGET_H

#include <iostream>
#include <iomanip>
#include <vector>


namespace single_photon
{
  //header printer; it returns the spacing of each columns;
  std::vector<int> Printer_header( std::vector< std::string> headings);

  //table printer
  void Printer_content( std::vector< std::string > nums, std::vector<int> spacers);

}

#endif
