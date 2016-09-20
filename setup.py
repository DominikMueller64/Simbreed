# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


#with open('README.rst') as f:
#    readme = f.read()
#
#with open('LICENSE') as f:
#    license = f.read()

setup(
    name='simbreed',
    version='0.1',
    description='Package for simulation in plant breeding.',
#    long_description=readme,
    author='Dominik Mueller',
    author_email='Dominik_Mueller@uni-hohenheim.de',
    url='https://github.com/X',
    license='GPLv3',
    packages=find_packages(exclude=('tests', 'docs', 'profiling', 'graveyard'))
)
