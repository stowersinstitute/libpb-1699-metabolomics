# Required Software

* GHC 8.6.5 (or later) and Cabal 2.4.1.0 (or later)
 - Can be obtained from [ghcup](https://www.haskell.org/ghcup/).
* Python 3.8.3 or later (if running 3.9 or later see below)
* pipenv 2018.11.26 or later
 - Can usually be obtained by running `sudo pip install pipenv` or `conda install pipenv` if using Anaconda.
* R (specifically the `Rscript` command) 3.6.3 or later with the bayesglm package installed

# Running

Run the following three commands:

```
cabal new-update
cabal new-build
cabal new-run pipeline
```

This will generate all figures under `out/fig`, and supplementary information under `out/supp`.

# Notes

## Newer Python Versions (3.9+)

If running a newer version of Python, you may need to edit the Pipfile and change the line:

```
python_version = "3.8"
```

to (substitute your Python version)

```
python_version = "3.9"
```

However, some packages may not work with newer versions of Python or may have different behavior.

## Outliers

Most analyses produce outputs with and without the following outliers:

* Pachon Muscle Refed 5
* Tinaja Liver Refed 6
* Pachon Liver 30d Starved 3

Results including the outliers are in an `outliers` subdirectory whereas results excluding the outliers are in a `no-outliers` subdirectory (this would be more clear if English had a word like *kein*, but alas it does not). Most results from the paper exclude outliers.
