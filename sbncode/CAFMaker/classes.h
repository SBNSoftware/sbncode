#include "sbnanaobj/StandardRecord/StandardRecord.h"

#include "lardataobj/RecoBase/Slice.h"
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// We have to use this syntax, introducing a dummy struct, instead of the
// simpler "template class" to avoid the need to define
// StandardRecord::operator<(). See https://cdcvs.fnal.gov/redmine/issues/7075
namespace{
  struct dict{
    art::PtrVector<caf::StandardRecord> m1;
  };
}
