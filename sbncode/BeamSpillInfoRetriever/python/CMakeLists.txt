# Analysis-level BNB beam-quality tools (pure Python, see README.md).
# Installed together, with their data, in <product>/python/bnb_beam_quality, since the
# scripts import each other and logger_fill.py reads acnet/ next to itself.

set(bq_scripts
  bnb_beam_quality.py
  dump_bnb_quality.py
  merge_bnb_dumps.py
  bnb_fom.py
  bnb_fom_fast.py
  neighbor_fill.py
  fom_fixups.py
  logger_fill.py
  validate_fom.py
  compare_cpp_caf.py
  spill_event_map.py
  plot_fom_compare.py
)

install(PROGRAMS ${bq_scripts} DESTINATION python/bnb_beam_quality)
install(FILES README.md requirements.txt DESTINATION python/bnb_beam_quality)
install(PROGRAMS acnet/d44_search.py DESTINATION python/bnb_beam_quality/acnet)
install(FILES acnet/README.md acnet/logger_fill_calib.json acnet/recovered_1min_backup_node96.csv
        DESTINATION python/bnb_beam_quality/acnet)
