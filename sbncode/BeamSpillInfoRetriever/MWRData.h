#ifndef _MWRDATA_H
#define _MWRDATA_H

#include <string.h>
namespace sbn{
class MWRData
{
  typedef struct swicheader_t {
    long timesec;
    long timensec;
    long gpstime1;
    long gpstime2;
    short boosterevent;
    short mievent;
    short hz15micnt;
    long delta1f;
    short pulsemi;
    short pulsesc;
  } swicheader_t;
  
  typedef struct mwrpulse_t {
    short hor[48];
    short ver[48];
    swicheader_t sheader;
  } mwrpulse_t;  
  
  static long flipByte(long data)
  {
    return ((data>>16)&0x0000FFFF) | ((data<<16)&0xFFFF0000);
  }
  
  mwrpulse_t getMWRdata(short* data, int nblock) const;
  
 public:
  std::vector< std::vector < int > > unpackMWR(std::string packed_data, std::vector<double> &time_stamp, double timeoffset=0) const;
};
}

#endif /* #ifndef _MWRDATA_H */
