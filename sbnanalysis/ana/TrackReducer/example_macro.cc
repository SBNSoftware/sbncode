#include "TFile.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "Data.h"

void read_file() {
  // open the input file
  TFile f("./data.root");

  // maximum number of events to read from file
  int max_event = -1;

  // setup a "TTreeReader" to read from the "sbnana" Tree
  TTreeReader myReader("sbnana", &f);

  // setup to read the branch "tracks"
  // Each instance on the branch has a
  // "Tracks" object which will have all the data on
  // the event that you need
  TTreeReaderValue<Tracks> tracks(myReader, "tracks");

  TFile out("output.root", "CREATE");

  // histogram which will store muon energies
  TH1F *muon_energy = new TH1F("muon_energy", "muon_energy", 50, 0, 5);

  // histogram which will store distance from reconstructed start
  // to true start
  TH1F *start_distance = new TH1F("start_distance", "start_distance", 50, 0, 5);

  int event_ind = 0;
  // Loop over the sbnana Tree
  while (myReader.Next()) {
    // break if we are done
    if (max_event >= 0 && event_ind == max_event) break;

    // We are going to make two plots:
    // A distribution of the true muon energies
    // and a distribution of the distance from the true 
    // muon start point to the reconstructed muon start point

    // The input simulation generates a single muon in the 
    // detector which then goes on to create a host of other
    // particles. To distinguish the initial particle
    // (which is the one we want) it is given a process name
    // of "primary"
    //
    // So to find our muon, we can search through the list
    // of particles for one with the process name "primary"

    int primary_index = -2;
    for (unsigned i = 0; i < tracks->truth.size(); i++) {
      const TrueParticle &particle = tracks->truth[i];
      if (particle.process == "primary") {
        std::cout << "energy: " << particle.energy << std::endl;
        muon_energy->Fill(particle.energy);    
        primary_index = i;
        break;
      }
    }

    for (const Track &track: tracks->reco_tracks) {
      // check that this is the track that matches
      // to the true primary muon
      if (track.truth_match == primary_index) {
        TVector3 true_start = &tracks->truth[primary_index].trajectory[0][0]; // convert array to TVector
        TVector3 track_start = &track.trajectory[0][0];
        start_distance->Fill((true_start - track_start).Mag());
        break;
      }
    }

    event_ind += 1;
  }
  out.cd();
  muon_energy->Write();
  start_distance->Write();
}
