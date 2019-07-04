#pragma once

#include <ctime>
#include <string>

namespace ana
{
  /// A simple ascii-art progress bar
  class Progress
  {
  public:
    /// Create and draw the progress bar
    Progress(const std::string& title);
    ~Progress();
    /// Update the progress fraction between zero and one
    void SetProgress(double frac);
    /// Call this when action is completed
    void Done();
  protected:
    std::string FormatTime(double sec) const;

    bool fDone; ///< Has \ref Done been called?
    int fIFrac; ///< What character are we on? Prevents unnecessary redraws

    time_t fStart;
    time_t fPrevCall;

    // Only one bar may be live at a time. Prevents overdrawing
    bool fLive; ///< Is this bar live (drawable?)
    static bool fAnyLive; ///< Are any bars live?
  };
}
