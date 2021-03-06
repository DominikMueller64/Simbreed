# -*- coding: utf-8 -*-

from setuptools import setup, find_packages


#with open('LICENSE') as f:
#    license = f.read()

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
    name='simbreed',
    version='0.1',
    description='Package for simulations in plant breeding.',
    long_description=readme(),
    author='Dominik Mueller',
    author_email='dominikmueller64@yahoo.de',
    url='https://github.com/DominikMueller64/Simbreed',
    license='GPLv3',
    packages=find_packages(exclude=('tests', 'docs', 'profiling')),
    install_requires=['uuid',
                      'shortuuid',
                      'sortedcollections',
                      'sortedcontainers'],
    include_package_data=True,
    zip_safe=False
)
