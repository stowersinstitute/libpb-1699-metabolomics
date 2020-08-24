# Required Software

* GHC 8.6.5 (or later) and Cabal 2.4.1.0 (or later)
 - Can be obtained from [ghcup](https://www.haskell.org/ghcup/).
* Python 3.8.3 or later (if running 3.9 or later see below)
* pipenv 2018.11.26 or later
 - Can usually be obtained by running `sudo pip install pipenv` or `conda install pipenv` if using Anaconda.
* R (`Rscript`) 3.6.3 or later with the bayesglm package installed

# Running

```
runhaskell shake.hs
```

This will generate all figures under `out/fig`.

# Newer Python Versions (3.9+)

If running a newer version of Python, you may need to edit the Pipfile and change the line:

```
python_version = "3.8"
```

to (substitute your Python version)

```
python_version = "3.9"
```
