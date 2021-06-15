from pathlib import Path


from setuptools import setup, find_packages

setup(
    name='cosg',
    version='1.0.0',
    description='Marker gene identification for single-cell sequencing data',
    long_description=Path('README.rst').read_text('utf-8'),
    url='https://github.com/genecell',
    download_url='https://github.com/genecell',
    packages=find_packages(exclude=['bin', 'conf', 'data', 'target']),
    install_requires=[
        l.strip() for l in Path('requirements.txt').read_text('utf-8').splitlines()
    ],
    include_package_data=True,
    author='Min Dai',
    author_email='daimin@zju.edu.cn',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Framework :: Jupyter',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ],
)

### Refer to https://github.com/theislab/scanpy/blob/master/setup.py