
import os
import setuptools
from distutils.core import setup

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='EXFOR_chart',
    version='0.0.1',
    author='Anatoli Fedynitch',
    author_email='afedynitch@gmail.com',
    py_modules=['exfor_chart.py', 'exfor_tools.py'],
    url='https://github.com/afedynitch/EXFOR_chart',
    install_requires=[
        'x4i3',
        'six',
        'periodictable',
        'numpy',
        'matplotlib'
    ],
    # extras_require={"tests": ["pytest"]},
    description='Interactive visulaization of photo-nuclear cross sections.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Physics',
        'Intended Audience :: Science/Research',
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: BSD License',
    ]
)
