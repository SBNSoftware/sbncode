#include "DynamicSelector.h"
#include "../uScript/api.h"
#include <iostream>

template <typename T>
void Increment(std::vector<unsigned> &indices, const std::vector<std::vector<T>> &lists) {
  for (unsigned i = 0; i < indices.size(); i++) {
    indices[i] += 1;
    if (indices[i] == lists[i].size()) {
      indices[i] = 0;
    }
    else {
      break;
    }
  }
}

bool NotZero(const std::vector<unsigned> &indices) {
  for (unsigned i = 0; i < indices.size(); i++) {
    if (indices[i] != 0) return true;
  }
  return false;
}

std::vector<std::string> numu::MultiplyNames(const std::vector<std::vector<std::string>> &strings, char join) {
  std::vector<std::string> ret;
  std::vector<unsigned> indices (strings.size(), 0);
  do {
    std::string this_str;
    for (unsigned i = 0; i < indices.size(); i++) {
      this_str += strings[i][indices[i]] + join;
    }
    ret.push_back(this_str);
    Increment(indices, strings);
  } while (NotZero(indices));
  return ret;
}

std::vector<numu::TrackSelector> numu::MultiplyTrackSelectors(const std::vector<std::vector<std::string>> &track_function_strings) {
  std::vector<numu::TrackSelector> ret;
  std::vector<std::vector<numu::TrackFunction>> track_functions;
  for (unsigned i = 0; i < track_function_strings.size(); i++) {
    track_functions.emplace_back();
    for (unsigned j = 0; j < track_function_strings[i].size(); j++) {
      track_functions[i].push_back(uscript::compile<numu::RecoTrack, numu::TrueParticle, unsigned>("track", "particle", "mctype", track_function_strings[i][j].c_str()));
    }
  }
  std::vector<unsigned> indices (track_functions.size(), 0);
  do {
    std::vector<numu::TrackFunction> functions;
    for (unsigned i = 0; i < indices.size(); i++) {
      functions.push_back(track_functions[i][indices[i]]);
    }
    ret.push_back(
      [functions](const numu::RecoTrack &track, const numu::TrueParticle &part, unsigned mctype) {
        for (const numu::TrackFunction &function: functions) {
          uscript::Value ret = function(&track, &part, &mctype);
          if (!ret) return false;
        }
        return true;
      });
    Increment(indices, track_functions);
  } while (NotZero(indices));
  return ret;
}
