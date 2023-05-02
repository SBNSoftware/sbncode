source /cvmfs/icarus.opensciencegrid.org/products/icarus/setup_icarus.sh
setup pythia8 v8_3_05 -q e20:p392:prof
python -m venv env
source env/bin/activate
pip install --upgrade pip
pip install numpy jupyterlab nbstripout tqdm ipywidgets
