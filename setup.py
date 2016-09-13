
from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='Scoary',
      version='1.5.0',
      description='Microbial pan-GWAS using the output from Roary',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GPLv3',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Medical Science Apps.',
        'Intended Audience :: Science/Research',
      ],
      keywords='microbial genomics GWAS Roary',
      url='https://github.com/AdmiralenOla/Scoary',
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
          'console_scripts': ['scoary=scoary.methods:main'],
      },
      include_package_data=True,
      zip_safe=False)
