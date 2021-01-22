////////////////////////////////////////////////////////////////////////
// \file    SRFlashMatch.h
// \brief   SRFlashMatch object for flashmatch summary information.
// \author  $Author: psihas@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef SRFLASHMATCH_H
#define SRFLASHMATCH_H


namespace caf
{
  /// A matching of TPC flashmatch charge to Optical flash light
  class SRFlashMatch
    {
    public:

      SRFlashMatch();
      virtual ~SRFlashMatch();

      bool  present;
      float score;
      float time;
      float pe;

      void setDefault();


    };
} // end namespace

#endif // SRFLASHMATCH_H
//////////////////////////////////////////////////////////////////////////////
