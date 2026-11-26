/**
 * @file   sbncode/CAFMaker/FillQELambda.h
 * @brief  CAFMaker utilities for filling "qelambda" objects into CAF.
 * @author Sergio Yebes (syebes@ugr.es)
 */

#ifndef FILLQELAMBDA_H
#define FILLQELAMBDA_H

#include "sbnobj/SBND/QELambdaReco/LambdaProducts.hh"
#include "sbnobj/SBND/QELambdaReco/LambdaPID.hh"
#include "sbnanaobj/StandardRecord/SRQELambda.h"

namespace caf
{

  void FillFRANSLambda( const qelambda::QELambdaReco LAr_LambdaRecoProducts, caf::SRQELambda& CAF_LambdaRecoProducts);
  void FillLambdaPID( const qelambda::QELambdaPID LAr_LambdaRecoProducts, caf::SRQELambda& CAF_LambdaRecoProducts);


}

#endif // FILLQELAMBDA_H