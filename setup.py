#!/usr/bin/env python

from distutils.core import setup

setup(
    name='lftc',
    version='1.0',
    description='Limit Flux To Core',
    author='The Quantitative Metabolic Modeling group',
    author_email='tbackman@lbl.gov',
    url='https://github.com/JBEI/LimitFluxToCore',
    packages=['lftc'],
    install_requires=['cobra'],
    license='see LICENSE file',
    download_url = 'https://github.com/JBEI/LimitFluxToCore/archive/1.0.tar.gz', 
    keywords = ['metabolism', 'flux'],
    classifiers = []
    )
