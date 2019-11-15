#!usr/bin/env python
from setuptools import setup

setup(name='APOExptime',
      version='0.0.1',
      author='Manuel H. Canas, Alexander Stone-Martinez, Bryson Stemock, Rogelio Ochoa, Hasan Rahman',
      author_email = 'canasmh@nmsu.edu',
      packages=['APOExptime'],
      scripts=['APOExptime/bin/etc_script.py'],
      package_data={'APOETC':['./data/apo3_5m/*.txt',
                              './data/apo3_5m/Arctic/*.dat',
                              './data/apo3_5m/Agile/*.dat',
                              './data/apo3_5,/Test/*.dat',
                              './data/Sky/*.txt']},
      description='An exposure time calculator for the ARC 3.5m telescope',
      include_package_data=True
      )
