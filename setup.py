#!/usr/bin/env python

from setuptools import setup

setup(name='limitsofgroth',
      entry_points = {
          'console_scripts' :
              ['limitsofgrowth = limitsofgrowth.main:main',
               ]},
      version='0.0.0',
      description='A simple World Model, after J.W. Forrester',
      author='Martin Ortbauer',
      author_email='mortbaue@gmail.com',
      url='http://github.com/mortbauer/limitsofgrowth/',
      download_url='http://github.com/mortbauer/limitsofgrowth/',
      license='MIT',
      packages=['limitsofgrowth'],
      install_requires=['pyqtgraph'],
      package_data={'limitsofgrowth':['data/default.conf']},
      classifiers=[
        'Development Status :: Alpha',
        'Topic :: Text Processing :: Markup',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Documentation',
        'License :: OSI Approved :: GNU General Public License (GPL)'],
)
