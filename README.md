# lftc
JBEI Limit Flux To Core Library 

# Installing
```
pip install -e .
```

# Build Sphinx documentation located in docs/_build/html/index.html
Please use Google style docstrings: http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
```
cd docs
make html
```

# Running test notebooks inside the Debian Cheminformatics Docker container
```
docker run -it --rm -v `pwd`:/f -w /f -p 8888:8888 tbackman/debian-cheminformatics jupyter notebook --no-browser --ip=* --allow-root
```
