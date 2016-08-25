#!/usr/bin/env python

from distutils.core import setup

setup(name='mutation-util',
    version='0.4.0',
    description='Python tools to identify somatic mutations.',
    author='Ken-ichi Chiba',
    author_email='kchiba@hgc.jp',
    url='https://github.com/Genomon-Project/GenomonFisher',
    package_dir = {'': 'lib'},
    packages=['mutil'],
    scripts=['mutil'],
    license='GPL-3'
)
