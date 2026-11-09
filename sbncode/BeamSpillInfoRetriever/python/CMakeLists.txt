# update_flatcafs.py and the BNB beam-quality code it uses (Python, see README.md).
# Installed together, with their data, in <product>/python/bnb_beam_quality, since the
# scripts import each other and logger_fill.py reads acnet/ next to itself.

set(bq_scripts
  update_flatcafs.py
  bunch_rotation.py
  dump_bnb_quality.py
  merge_bnb_dumps.py
  bnb_fom.py
  bnb_fom_fast.py
  neighbor_fill.py
  fom_fixups.py
  logger_fill.py
)

install(PROGRAMS ${bq_scripts} DESTINATION python/bnb_beam_quality)
install(FILES README.md requirements.txt DESTINATION python/bnb_beam_quality)
install(FILES acnet/README.md acnet/logger_fill_calib.json acnet/recovered_1min_backup_node96.csv
        DESTINATION python/bnb_beam_quality/acnet)
