
from setuptools import setup
from scoary import __version__ as sv


def readme():
    with open('README_pypi.md') as f:
        return f.read()


setup(name='scoary',
      version=sv,
      description='Microbial pan-GWAS using the output from Roary',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Intended Audience :: Science/Research',
      ],
      keywords=['microbial', 'genomics', 'GWAS', 'Roary'],
      url='https://github.com/AdmiralenOla/Scoary',
      download_url='https://github.com/AdmiralenOla/Scoary/tarball/v1.5.0',
      author='Ola Brynildsrud',
      author_email='ola.brynildsrud@fhi.no',
      license='GPLv3',
      packages=['scoary'],
      install_requires=[
          'scipy>=0.16',
          'argparse'
      ],
      test_suite='nose.collector',
      tests_require=[],
      entry_points={
          'console_scripts': ['scoary=scoary.methods:main',
                             'scoary_GUI=scoary.GUI:main'],
      
      },
      include_package_data=True,
      zip_safe=False)
