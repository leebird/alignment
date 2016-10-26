#!/usr/bin/env python

import os
import setuptools

VERSION = '0.0.1'

README = open(os.path.join(os.path.dirname(__file__), 'README.md')).read()

setuptools.setup(
    name='alignment',
    author='Gang Li',
    author_email='leemagpie@gmail.com',
    url='https://github.com/leebird/alignment',
    description='Alignment algorithms for aligning text.',
    long_description=README,
    py_modules=['alignment'],
    license='BSD',
    # test_suite='tests',
    version=VERSION,
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 1 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ],
    platforms='any',
)
