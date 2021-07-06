#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include "MWRData.h"

using namespace std;

namespace sbn{

  std::vector< std::vector < int > > MWRData::unpackMWR(std::string packed_data, std::vector<double> &time_stamp, double timeoffset) const
{

  std::vector<std::vector<int> > unpacked_data;
  unpacked_data.resize(4);
  short data[444];

  std::vector<std::string> row(0);
  boost::split(row, packed_data, boost::is_any_of(","));
  if (row.size()==447) {
    for (int i=3;i<447;i++) {
      data[i-3]=atoi(row[i].c_str());
    }
    string devname=row[1].substr(0,8);
    for (int idev=0;idev<4;idev++) {
      mwrpulse_t mwr=getMWRdata(data,idev);
      time_stamp.push_back(mwr.sheader.timesec+mwr.sheader.timensec/1000000000.+timeoffset);
      for (int ich=0;ich<48;ich++) {
	unpacked_data[idev].push_back(mwr.hor[ich]);
      }
      for (int ich=0;ich<48;ich++) {
	unpacked_data[idev].push_back(mwr.ver[ich]);
      }
    }
  } else {
    cout <<"Bad data!"<<endl;
    return unpacked_data;
  }

  return unpacked_data;
}

MWRData::mwrpulse_t MWRData::getMWRdata(short* data, int nblock) const
{
  mwrpulse_t mwrdata;

  memcpy(&mwrdata.hor,                  data+nblock*111,    96);
  memcpy(&mwrdata.ver,                  data+nblock*111+48, 96);
  memcpy(&mwrdata.sheader.timesec,      data+nblock*111+96,  4);
  memcpy(&mwrdata.sheader.timensec,     data+nblock*111+98,  4);
  memcpy(&mwrdata.sheader.gpstime1,     data+nblock*111+100, 4);
  memcpy(&mwrdata.sheader.gpstime2,     data+nblock*111+102, 4);
  memcpy(&mwrdata.sheader.boosterevent, data+nblock*111+104, 2);
  memcpy(&mwrdata.sheader.mievent,      data+nblock*111+105, 2);
  memcpy(&mwrdata.sheader.hz15micnt,    data+nblock*111+106, 2);
  memcpy(&mwrdata.sheader.delta1f,      data+nblock*111+107, 4);
  memcpy(&mwrdata.sheader.pulsemi,      data+nblock*111+109, 2);
  memcpy(&mwrdata.sheader.pulsesc,      data+nblock*111+110, 2);

  mwrdata.sheader.timesec=flipByte(mwrdata.sheader.timesec);
  mwrdata.sheader.timensec=flipByte(mwrdata.sheader.timensec);
  mwrdata.sheader.gpstime1=flipByte(mwrdata.sheader.gpstime1);
  mwrdata.sheader.gpstime2=flipByte(mwrdata.sheader.gpstime2);
  mwrdata.sheader.delta1f=flipByte(mwrdata.sheader.delta1f);

  return mwrdata;

}
}
