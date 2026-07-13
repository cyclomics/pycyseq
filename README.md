<p>
  <img src="assets/header.png" alt="header"/>
</p>

# pycyseq

`pycyseq` is a python module to read and manipulate CySeq data.

For documentation please visit: 

## Installation

We recommend installing `pycyseq` in a python virtual environment. This can be easily done in python itself, or by using a virtual environment manager such as [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/getting-started.html#managing-python) or [uv](https://docs.astral.sh/uv/).


```bash
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

## Development

For developers, please install [uv](https://docs.astral.sh/uv/getting-started/installation/).

```bash
git clone https://github.com/cyclomics/pycyseq
cd pycyseq
uv sync --dev
```

### Useful commands

**Run tests:** `uv run pytest`

**Preview documentation:** `uv run zensical serve`

**Check typing:** `uv run ty check`

**Check formatting:** `uv run ruff format && uv run ruff check`

### License and Copyright

(c) 2026 Cyclomics

pycyseq is is distributed under the terms of the Cyclomics Research Source License, v. 1.0. If a copy of the License was not distributed with this file, You can obtain one at https://github.com/cyclomics/pycyseq.