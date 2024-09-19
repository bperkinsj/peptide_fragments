A series of scripts to match domains and secondary structures to mass spectrometry data.

## Cloning the Repository

You can follow the instructions [here](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) for cloning a repository.

## Installation

First make sure that you have [Anaconda installed](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Navigate to this repo in the terminal and [install the virtual environment](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html):

```
$ conda env create -f environment.yml
$ conda activate peptides
```

## Input data
Place your files of interest into the 'input' folder (or whatever directory you have assigned to 'input' in config.ini).

Make sure that you don't have multiple values in any columns such as '20,40' in PEP.PeptidePosition or 'P28482/P60428' in PG.ProteinAccessions. The script 'tidy_data.py' will fix columns with values separated by ',' or ';' -  just pass it the directory containing your files:

```
$ python tidy_data.py --dir='input'
```

## Running

Activate the environment

```
$ conda activate peptides
```

To run, call:

```
$ python main.py --config configs/config.ini
```

Config.ini contains a few variables you can change. In addition, if you would like to run only a single file, you can either adjust it in config.ini or pass it to main.py via the --input command (note that you'll also have to pass the --is_file flag either way):

```
$ python main.py --config configs/config.ini --input='input/file.txt' --is_file
```

The script automatically deletes any downloaded AlphaFold2 files at its completion, but if you would like to keep them you can pass the --keep_structures flag when you call the script.

## Results

Results can be found in 'data/[experiment_name]/results. 

The "alpha_helix", "beta_sheet", and "bend" columns show residues that are noted to have these features in the structure files. The "type" column notes what kind of domain feature is present. The "region" column shows the region demarcating the domain feature. The "description" column gives a more detailed description of what the domain feature is.