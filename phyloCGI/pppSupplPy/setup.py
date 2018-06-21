from setuptools import setup, find_packages

setup(
    name = 'pppSupplPy',
    version = '0.0.2',
    py_modules = ['pppSupplPy'],

    author = 'YulongNiu',
    author_email = 'yulong.niu@hotmail.com',

    url = 'http://yulongniu.bionutshell.org/',
    description = 'A supplement Python package for PrePhyloPro',

    packages = find_packages(),
    
     project_urls={
        'Web application': 'http://bioinfor.scu.edu.cn/phyloprofile/',
        'Source': 'https://github.com/YulongNiu/PhyloProfileWeb',
    },
)

