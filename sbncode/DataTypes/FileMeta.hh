#ifndef __sbnanalysis_core_FileMeta__
#define __sbnanalysis_core_FileMeta__

/**
 * \file FileMeta.hh
 *
 * File-level information.
 *
 * Author: G. Putnam <grayputnam@uchicago.edu> 2019/12/11
 */

/**
 * \class FileMeta
 * \brief Metadata for each input file.
 */
class FileMeta {
public:
  FileMeta()
      : n_events(0), n_gen_events(0) {}

  FileMeta(unsigned _n_events, unsigned _n_gen_events)
    : n_events(_n_events), n_gen_events(_n_gen_events) {}

  unsigned n_events;
  unsigned n_gen_events;
};

#endif  // __sbnanalysis_core_FileMeta__

