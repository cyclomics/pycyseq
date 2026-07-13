# Installation

We recommend installing `pycyseq` in a python virtual environment. This can be easily done in python itself, or by using a virtual environment manager such as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-python) or [uv](https://docs.astral.sh/uv/).


```bash title="Example pycyseq installation in a python virtual environment"
# download the pycyseq repository and move in it
git clone https://github.com/cyclomics/pycyseq
cd pycyseq

# create the python virtual environment, and activate the environment
python3 -m venv venv
source venv/bin/activate

# update pip dependencies and install pycyseq
python3 -m pip install --upgrade pip
pip3 install . --no-cache-dir
```