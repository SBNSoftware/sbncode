#include "CAFAna/Core/Progress.h"

#include "sys/stat.h"

#include "TString.h"

#include <iostream>

namespace{
  const int kBarWidth = 60;
}

namespace ana
{
  bool Progress::fAnyLive = false;

  //----------------------------------------------------------------------
  Progress::Progress(const std::string& title)
    : fDone(false), fIFrac(-1), fStart(time(0)), fPrevCall(time(0)), fLive(false)
  {
    // If no else is drawing, we can
    if(!fAnyLive){
      fLive = true;
      fAnyLive = false;
    }

    if(!fLive) return;

    std::cout << title << "..." << std::endl;
    SetProgress(0); // Draw the initial bar
  }

  //----------------------------------------------------------------------
  Progress::~Progress()
  {
    // Finish up in case the user forgot to call Done()
    Done();

    // If we were the ones drawing, we're not anymore
    if(fLive) fAnyLive = false;
  }

  //----------------------------------------------------------------------
  void Progress::SetProgress(double frac)
  {
    if(!fLive || fDone) return;

    // Check if we're outputting to a file. If so don't bother showing off
    // with the progress bar, it won't work.
    struct stat buf;
    fstat(fileno(stdout), &buf);
    const bool isFile = (buf.st_mode & S_IFREG) || (buf.st_mode & S_IFIFO);
    if(isFile) return;

    const int ifrac = (kBarWidth-1)*frac;

    const time_t t_now = time(0);

    // Don't repaint unnecessarily
    if(ifrac == fIFrac && t_now - fPrevCall < 2) return;

    fIFrac = ifrac;
    fPrevCall = time(0);

    std::string str(kBarWidth, ' ');
    for(int i = 0; i < ifrac; ++i) str[i] = '=';
    str[ifrac] = '>';
    str[0] = '[';
    str[kBarWidth-1] = ']';

    if(frac > 0){
      const int elapse = t_now - fStart;
      if(elapse > 2){ // Don't show for very short steps
        if(frac < 1)
          str += " "+FormatTime(elapse*(1-frac)/frac);
        else
          str += " "+FormatTime(elapse);
        str += "    "; // Enough to cover up any previous version
      }
    }

    std::cout << "\r" << str << std::flush;

    if(frac == 1){
      fDone = true;
      std::cout << std::endl;
    }
  }

  //----------------------------------------------------------------------
  void Progress::Done()
  {
    if(!fLive) return;

    if(fDone) return; // Can easily be called multiple times

    SetProgress(1); // Make sure the bar shows 100%

    fDone = true;
  }

  //----------------------------------------------------------------------
  std::string Progress::FormatTime(double sec) const
  {
    // Yes, I'm sure there's a standard way to do this, but this was easy, and
    // lets me print exactly what I want.
    std::string ret;
    if(sec >= 60*60-.5){
      ret += TString::Format("%dh", (int(sec+.5)/(60*60))).Data();
    }
    if(sec >= 60-.5){
      ret += TString::Format("%dm", (int(sec+.5)/60)%60).Data();
    }
    if(sec < 60*60){ // don't clutter if still measured in hours
      ret += TString::Format("%ds", (int(sec+.5)%60)).Data();
    }
    return ret;
  }
}
