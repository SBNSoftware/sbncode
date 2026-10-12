# ACNet logger readings of BNB BPMs missing from IFBeam (ICARUS Run 2, 2023)

Used by `logger_fill.py` (`update_flatcafs.py --logger-fill`). The D44 query script that recovered them is in the code archive of the beam-quality technote.


Result: the per-spill (15 Hz event) data is not kept on any node. 1-minute readings for HP875, VP875, HPTG1, HPTG2 and VPTG2
do exist on node 96 (Backup), and node 70 (MBoon) returns the same values at p,60000. VP873 is not logged on any periodic node.

Recovery query, which gives every reading in the time range (times printed are US Central):
https://www-bd.fnal.gov/D44/dataLogger?t0=<ms>&t1=<ms>&how2plot=all&how_much=all&x=Time&y0=<DEV>&y0_node=96&y0_color=black&Y0_min=auto&Y0_max=auto&Y0_event=p,300000,true&txt=true
  Outage A (HP875/VP875): t0=1686344400000 t1=1686592800000
  Outage B (HPTG1/HPTG2/VPTG2): t0=1678906800000 t1=1678975200000

`recovered_1min_backup_node96.csv` holds the readings recovered from node 96 (Backup), one per minute, for both full outages. Columns: device, outage, utc_ms, utc_iso, value_mm.
Each UTC time is t0 from the query plus the offset in seconds that D44 reports. There are 3807 rows each for HP875 and VP875, and 1085 rows each for HPTG1, HPTG2 and VPTG2. VP873 has no rows because it isn't logged on Backup.
