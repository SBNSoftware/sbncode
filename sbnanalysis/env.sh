# Set up sbnanalysis

export SBN_ANALYSIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export LD_LIBRARY_PATH=$SBN_ANALYSIS_DIR/lib:$LD_LIBRARY_PATH

