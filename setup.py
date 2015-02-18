from setuptools import setup
from distutils.core import Extension
from Cython.Distutils import build_ext
import numpy as np

ext_dtram = Extension(
        "pytram.dtram.ext",
        sources=["ext/dtram/dtram.pyx", "ext/dtram/_dtram.c" ],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3"]
    )
ext_xtram = Extension(
        "pytram.xtram.ext",
        sources=["ext/xtram/xtram.pyx", "ext/xtram/_xtram.c" ],
        include_dirs=[np.get_include()],
        extra_compile_args=["-O3"]
    )

setup(
    cmdclass={'build_ext': build_ext},
    ext_modules=[
            ext_dtram,
            ext_xtram
        ],
    name='pytram',
    version='0.1.6',
    description='The TRAM package',
    long_description='The python interface to the TRAM framework for estimating Markov state models from biased MD simulations.',
    classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: C',
            'Programming Language :: Cython',
            'Programming Language :: Python :: 2.7',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics'
        ],
    keywords=[
            'free energy',
            'Markov state model',
            'TRAM',
            'dTRAM',
            'xTRAM'
        ],
    url='https://github.com/markovmodel/pytram',
    author='The pytram team',
    author_email='pytram@lists.fu-berlin.de',
    license='Simplified BSD License',
    setup_requires=[
            'numpy>=1.7.1',
            'cython>=0.15',
            'setuptools>=0.6'
        ],
    tests_require=[ 'numpy>=1.7.1', 'nose>=1.3' ],
    install_requires=[ 'numpy>=1.7.1' ],
    packages=[
            'pytram',
            'pytram.dtram',
            'pytram.xtram'
        ],
    test_suite='nose.collector',
    scripts=[
            'bin/dtram.py',
            'bin/xtram.py'
        ]
)
