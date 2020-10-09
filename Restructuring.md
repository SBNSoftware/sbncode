Repository restructuring
=========================


This is a temporary running document to summarize the changes happening in
SBN code and the related modifications.

The repositories currently involved are: `sbncode`, `sbnobj`, `sbndcode` and `icaruscode`.
The latter two become dependent on `sbncode`, which depends on `sbnobj`.


Code transfer
--------------

| from                          | to                                 | to commits  | missing steps | partial? |
| ----------------------------- | ---------------------------------- | ----------- | ------------- | -------- |
| `icaruscode/PMT/Trigger/Data` | `sbnobj/ICARUS/PMT/Trigger/Data`   | `0834ea71`  |               | _no_     |
| `icaruscode/CRT/CRTProducts`  | `sbnobj/Common/CRT`                | `629ef134`  |               | `CRTHit.hh`, `CRTHit.cc`, `CRTTrack.hh`, `CRTTrack.cc`, `CRTTzero.hh`, `CRTTzero.cc`, `classes_def.xml`, `classes.h` |
| `icaruscode/CRT/CRTProducts`  | `sbnobj/ICARUS/CRT`                | `40ba93a6`  |               | `CRTData.cc`, `CRTData.hh`, `classes_def.xml`, `classes.h` |
| `sbndcode/CRT/CRTProducts`    | `sbnobj/SBND/CRT`                  | `40ba93a6`  |               | `CRTData.cc`, `CRTData.hh`, `classes_def.xml`, `classes.h` |
| `icaruscode/IcarusObj`        | `sbnobj/Common/Analysis`           | `743e2f50`  |               | _no_     |


In addition, "new" files were created:
* in `sbnobj/Common/CRT`: `CRTHit_Legacy.hh`, `CRTTzero_Legacy.hh`, `CRTTrack_Legacy.hh` (see [SBNSoftware/sbnobj pull request #2](https://github.com/SBNSoftware/sbnobj/pull/2))



Things to do when moving code
------------------------------

1.  copy the code to the new location _[A]_
2.  rework the destination path if needed
3.  update the `CMakeLists.txt` file of the destination directory to include the new paths, if any _[B]_
4.  update the fixing script `sbncode/scripts/updates/restructuring.py` or create a new one
5.  add and commit into the feature branch
6.  remove the code to the old location (and commit into the feature branch) _[A]_
7.  purge the `CMakeLists.txt` file of the old directory from the old path, if needed _[B]_
8.  run the update script and commit again (separating the addition and the fix into distinct commits)
9.  check that everything compiles and automatic tests are ok
10. add a line to the table above

(_[A]_ and _[B]_ mark mirror steps)
