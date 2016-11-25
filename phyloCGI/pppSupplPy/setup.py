from distutils.core import setup
from setuptools import setup, find_packages

setup(
    name = 'pppSupplPy',
    version = '0.0.1',
    py_modules = ['pppSupplPy'],

    author = 'YulongNiu',
    author_email = 'niuylscu@gmail.com',

    url = 'http://yulongniu.bionutshell.org/',
    description = 'A supplement python package for PrePhyloPro',

    packages = find_packages(),
)

