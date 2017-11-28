#!/usr/bin/env python

from distutils.core import setup

setup(
    name='lftc',
    version='1.0',
    description='Limit Flux To Core',
    author='The Quantitative Metabolic Modeling group',
    author_email='tbackman@lbl.gov',
    url='https://github.com/JBEI/limitfluxtocore',
    packages=['lftc'],
    install_requires=['cobra', 'numpy'],
    license='see license.txt file',
    download_url = 'https://github.com/JBEI/limitfluxtocore/archive/1.0.tar.gz', 
    keywords = ['metabolism', 'flux'],
    classifiers = [],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    )
