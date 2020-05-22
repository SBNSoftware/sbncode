#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

void read_numu_file() {
  // open the input file
  TFile f("/sbnd/data/users/gputnam/NuMu/outputs/sbnd_passthru_v1/SBND_numu_all.root");

  // histogram which will store neutrino energies
  TH1F *hist = new TH1F("", "", 50, 0, 5);

  // maximum number of events to read from file
  int max_event = -1;

  // setup a "TTreeReader" to read from the "sbnana" Tree
  // Each entry in the branch is a single "Event" in SBND
  TTreeReader myReader("sbnana", &f);

  // Access to energy of each neutrino in the event.
  // In SBND, there are sometimes multiple neutrinos interacting
  // in the detector in each event, so we need to keep track of an
  // array of energies. However, most of the time this array will
  // just have one entry.
  TTreeReaderArray<double> neutrino_energies(myReader, "truth.neutrino.energy");

  int event_ind = 0;
  // Loop over the sbnana Tree
  while (myReader.Next()) {
    // break if we are done
    if (max_event >= 0 && event_ind == max_event) break;
 
    // fill each neutrino enrgy in each interaction
    for (double e: neutrino_energies) {
      hist->Fill(e);
    }

    event_ind += 1;
  }
  hist->Draw();

}
