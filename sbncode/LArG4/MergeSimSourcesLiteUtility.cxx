#include "MergeSimSourcesLiteUtility.h"

sbn::MergeSimSourcesLiteUtility::MergeSimSourcesLiteUtility(const std::vector<int>& offsets)
  : sim::MergeSimSourcesUtility(offsets)
{
  fG4TrackIDOffsets = offsets;
}

void sbn::MergeSimSourcesLiteUtility::MergeMCParticleLites(
    std::vector<sim::MCParticleLite>& merged_vector,
    const std::vector<sim::MCParticleLite>& input_vector,
    size_t source_index)
{
  if (source_index >= fG4TrackIDOffsets.size())
    std::runtime_error("ERROR in MergeSimSourcesLite: Source index out of range!");

  //fMCParticleLiteListMap[source_index].resize(input_vector.size());
  merged_vector.reserve(merged_vector.size() + input_vector.size());

  std::pair<int, int> range_trackID(std::numeric_limits<int>::max(),
                                    std::numeric_limits<int>::min());

  for (size_t i_p = 0; i_p < input_vector.size(); ++i_p) {
    // sim::MCParticleLite does not have (yet?) the nice copy constructor with offset
    // like simb::MCParticle does. Hence we must do it the pedestrian way.
    //merged_vector.emplace_back(input_vector[i_p], fG4TrackIDOffsets[source_index]);
    merged_vector.emplace_back(input_vector[i_p]);
    auto offset = fG4TrackIDOffsets[source_index];
    auto tid = input_vector[i_p].TrackID() >= 0 ? input_vector[i_p].TrackID() + offset : input_vector[i_p].TrackID() - offset;
    merged_vector.back().TrackID(tid);

    //fMCParticleLiteListMap[source_index][i_p] = merged_vector.size() - 1;

    //if (std::abs(merged_vector.back().TrackID()) < range_trackID.first)
    //  range_trackID.first = std::abs(merged_vector.back().TrackID());
    //if (std::abs(merged_vector.back().TrackID()) > range_trackID.second)
    //  range_trackID.second = std::abs(merged_vector.back().TrackID());

    //UpdateG4TrackIDRange(range_trackID, source_index);
  }
}
