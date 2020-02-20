#ifndef CAF_CAFMAKERPARAMS_H
#define CAF_CAFMAKERPARAMS_H

#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/OptionalTable.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "canvas/Utilities/InputTag.h"
#include "art/Framework/Core/EDAnalyzer.h"

namespace caf
{
  struct CAFMakerParams
  {
    template<class T> using Atom = fhicl::Atom<T>;
    template<class T> using Table = fhicl::Table<T>;
    using Comment  = fhicl::Comment;
    using Name     = fhicl::Name;
    using string   = std::string;
    using InputTag = art::InputTag;

    /* Atom<bool> EnableBlindness */
    /* { */
    /*   Name("EnableBlindness"), */
    /*   Comment("true = hide sensitive info, false = include full record") */
    /* }; */

    Atom<std::string> CAFFilename { Name("CAFFilename"),
      Comment("Provide a string to override the automatic filename.")
    };

    Atom<string> ClusterLabel    { Name("ClusterInput") };
    Atom<string> DataTier        { Name("DataTier") };
    Atom<string> FileExtension   { Name("FileExtension") };
    Atom<string> GeneratorLabel  { Name("GeneratorInput") };

    Atom<bool> StrictMode        { Name("StrictMode"),
	Comment("Abort if any required product not found, unless label is empty")
    };

    Atom<bool> CutClearCosmic {
      Name("CutClearCosmic"),
      Comment("Cut slices which are marked as a 'clear-cosmic' by pandora"),
      true
    };

    Atom<bool> SelectOneSlice {
      Name("SelectOneSlice"),
      Comment("Only select one slice per spill (ranked by nu_score) [TODO: implement]."),
      false
    };

    fhicl::OptionalSequence<std::string> PandoraTagSuffixes {
      Name("PandoraTagSuffixes"),
      Comment("List of suffixes to add to TPC reco tag names (e.g. cryo0 cryo1)")
    };


  };
}

#endif
