#ifndef __sbnanalysis_core_Experiment__
#define __sbnanalysis_core_Experiment__

/**
 * \file Experiment.hh
 *
 * Identifiers for experiments.
 *
 * Author: A. Mastbaum <mastbaum@uchicago.edu>, 2019/03/01
 */

/** Identifier for known experiments. */
typedef enum {
  kExpSBND,
  kExpMicroBooNE,
  kExpICARUS,
  kExpDUNEND,
  kExpDUNEFD,
  kExpLArIAT,
  kExpOther = 1000
} Experiment;

#endif  // __sbnanalysis_core_Experiment__

