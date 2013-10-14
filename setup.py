#!/usr/bin/env python

from setuptools import setup
from setuptools.extension import Extension
from Cython.Distutils import build_ext

setup(name='limitsofgrowth',
      entry_points = {
          'console_scripts' :[
              'limitsofgrowth-qt = limitsofgrowth.qt_gui:main',
              'limitsofgrowth-kivy = limitsofgrowth.kivy_gui:main',
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
      package_data={'limitsofgrowth':
                    ['data/locale/de/LC_MESSAGES/limitsofgrowth.mo',
                     'data/sensitivityanalysis.png',
                     'data/world_indicators.hdf',
                     ]},
      ext_modules=[Extension('limitsofgrowth.traversedata',
                             ['limitsofgrowth/traversedata.pyx'])
                   ],
      cmdclass={"build_ext": build_ext},
      zip_safe=False,
      classifiers=[
        'Development Status :: Alpha',
        'Topic :: Text Processing :: Markup',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Topic :: Software Development :: Documentation',
        'License :: OSI Approved :: GNU General Public License (GPL)'],
)
