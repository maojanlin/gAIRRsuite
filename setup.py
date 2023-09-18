from distutils.core import setup
setup(
  name = 'gAIRR_suite',
  packages = ['gAIRR_suite', 'gAIRR_suite.scripts', 'gAIRR_suite.scripts.analysis', 'gAIRR_suite.material'],
  version = '0.2.0',
  license='MIT',
  description = 'Profiling genes encoding the adaptive immune receptor repertoire with gAIRR Suite.',
  author = 'Mao-Jan Lin',
  author_email = 'mj.maojanlin@gmail.com',
  url = 'https://github.com/maojanlin/gAIRRsuite',
  download_url = 'https://github.com/maojanlin/gAIRRsuite/tarball/master',
  keywords = ['AIRR', 'annotation', 'allele typing', 'novel allele'],
  install_requires=[
          'numpy',
          'pysam',
          'pyfastx'
      ],
  include_package_data=True,
  zip_safe = False,
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
  ],
  entry_points={"console_scripts": ["gAIRR_annotate = gAIRR_suite.gAIRR_annotate:main","gAIRR_call = gAIRR_suite.gAIRR_call:main"],},
)



