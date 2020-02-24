from distutils.core import setup
setup(
  name = 'APOExptime',         # How you named your package folder (MyLib)
  packages = ['src'],   # Chose the same as "name"
  version = '0.1',      # Start with a small number and increase it with every change you make
  license='AGPLv3',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Exposure time calculator for Apache point observatory',   # Give a short description about your library
  author = 'Alexander Stone-Martinez',                   # Type in your name
  author_email = 'stonemaa@nmsu.edu',      # Type in your E-Mail
  url = 'https://github.com/APOExposureTimeCalculator/APOExptime',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/APOExposureTimeCalculator/APOExptime/archive/v_01.tar.gz',    # I explain this later on
  keywords = ['Exposure time', 'Astronomy', 'APO'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'astropy',
          'yaml',
          'synphot',
          'scipy',
          'numpy'
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
    'Intended Audience :: Astronomers',      # Define that your audience are developers
    'Topic :: Astronomy Tools :: Exposure time',
    'License :: OSI Approved :: AGPLv3 License',   # Again, pick a license
    'Programming Language :: Python :: 3.6',
  ],
  include_package_data=True
)