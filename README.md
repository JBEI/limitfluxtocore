# JBEI Limit Flux To Core (lftc) Library 

This library provies a set of algorithms which 1) systematically calculate flux bounds for any specified "core" of a genome-scale model so as to satisfy the bow tie approximation and 2) automatically identify an updated and optimal set of core reactions that can satisfy this approximation more efficiently. First, we leverage linear programming to simultaneously identify the lowest fluxes from peripheral metabolism into core metabolism compatible with the observed growth rate and extracellular metabolite exchange fluxes. Second, we use Simulated Annealing to identify an updated set of core reactions that allow for a minimum of fluxes into core metabolism to satisfy these experimental constraints.

Please see license.txt and legal.txt for license details.

# Citing

If you use this software for published research, please cite the following journal article:
Backman, T.W.H.; Ando, D.; Singh, J.; Keasling, J.D.; García Martín, H. Constraining Genome-Scale Models to Represent the Bow Tie Structure of Metabolism for 13C Metabolic Flux Analysis. Metabolites 2018, 8, 3.  
http://www.mdpi.com/2218-1989/8/1/3

# Demo notebook

A demonstration usage notebook is provided at https://github.com/JBEI/limitfluxtocore/blob/master/notebooks/lftc_example.ipynb

# Installing
This software requires Python 3.x. You can install with pip after cloning and entering the git repo
with the following command. This may be called 'pip3' for the Python 3.x version on many systems.
```
pip install .
```

# Python API and documentation

See http://htmlpreview.github.io/?https://github.com/JBEI/limitfluxtocore/blob/master/docs/_build/html/index.html

# Use with jQMM

To use this library with the jQMM flux modeling library, first process and export your
model as shown in the above notebook. You must export the model with COBRApy as sbml2 as follows, and must have libSBML installed:
```
cobra.io.sbml.write_cobra_model_to_sbml_file(model, 'outputFile.xml', sbml_level=2, sbml_version=1,use_fbc_package=False)
```

Afterwards, this model can be used in jQMM by skipping the built in Limit Flux to Core step with the
limitFlux2Core=False option to findFluxesRanges(). You should also skip adding any new extracellular
exchange fluxes, as those should be added before running lftc.

# Instructions for lftc developers

## Build Sphinx documentation located in docs/_build/html/index.html
Please use Google style docstrings: http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
```
cd docs
make html
```

## Running test notebooks inside the Debian Cheminformatics Docker container
```
docker run -it --rm -v `pwd`:/f -w /f -p 8888:8888 tbackman/debian-cheminformatics jupyter notebook --no-browser --ip=* --allow-root
```
