# Set up sbnanalysis

export SBN_ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SBN_LIB_DIR=$SBN_ANALYSIS_DIR/build/lib
export LD_LIBRARY_PATH=$SBN_LIB_DIR:$LD_LIBRARY_PATH

