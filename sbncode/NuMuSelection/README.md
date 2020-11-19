This directory contains jupyter notebooks for processing _FLAT_ CAF files (.flat.root) 
in the nu-mu event selection. To set things up, make a python virtual
environment:

```
virtualenv env
. env/bin/activate
pip install -r requirements.txt
```

If `virtualenv` doesn't work, you can try running `~gputnam/virtualenv-15.1.0/virtualenv.py`

Then find a `.flat.root` CAF file and run the notebooks!

_NOTE_: before comitting anything to github, please strip the output of all jupyter notebooks.
You can do this by running:

```
nbstripout *.ipynb
```
