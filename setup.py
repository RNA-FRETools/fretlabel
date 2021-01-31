import setuptools

with open('README.md', 'r') as fh:
    long_description = fh.read()

about = {}
with open('fluordynamics/__about__.py') as a:
    exec(a.read(), about)

INSTALL_REQUIRES = [
    'numpy',
    'tqdm',
    'pandas',
    'biopandas',
    'pybind11',
    'jsonschema'
    ]

class get_pybind_include(object):
    def __str__(self):
        import pybind11
        return pybind11.get_include()

ext_modules = [setuptools.Extension('relaxation', sources=['fluordynamics/relaxation.cpp'], include_dirs=[get_pybind_include()], language='c++')]

setuptools.setup(
    name=about['__title__'],
    version=about['__version__'],
    author=about['__author__'],
    author_email=about['__email__'],
    description=about['__description__'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    url=about['__url__'],
    packages=setuptools.find_packages('pybind11', exclude=['fluordynamics/fluorlabel/gui.py']),
    install_requires=INSTALL_REQUIRES,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    keywords=about['__keywords__'],
    ext_modules=ext_modules,
    skripts=['gromacs-tools/solvate.sh', 'gromacs-tools/single_run.sh', 'gromacs-tools/continue_run.sh', 'gromacs-tools/resp_fit.sh'],
)
