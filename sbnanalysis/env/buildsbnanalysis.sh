
# Build sbnanalysis
cd $MRB_SOURCE/sbncode/sbnanalysis
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug
make install -j4
source bin/setup_sbnanalysis.sh

