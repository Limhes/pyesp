#!/usr/bin/env python3
# encoding: utf-8

#from distutils.core import setup, Extension
from setuptools import setup
from setuptools.extension import Extension

setup(
    name='PyESP',
    version='0.1.0',
    description='ESP24 port to Python 3',
    author='René Becker',
    author_email='limhes@gmail.com',
    url='http://limhes.net/PyESP',
    ext_modules=[
        Extension(
            'PyESP',
            sources = ['PyESP.c', 'simulate.c']
        )
    ]
)
