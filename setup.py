import setuptools
import pybind11
print(pybind11.get_include())

with open('README.md', 'r') as fh:
    long_description = fh.read()

about = {}
with open('fluordynamics/__about__.py') as a:
    exec(a.read(), about)

INSTALL_REQUIRES = [
    'numpy',
    'tqdm',
    'pybind11',
    'pandas',
    'biopandas'
    ]

relaxation = setuptools.Extension('relaxation', sources=['fluordynamics/relaxation.cpp'], include_dirs=[pybind11.get_include()], language='c++')

setuptools.setup(
    name=about['__title__'],
    version=about['__version__'],
    author=about['__author__'],
    author_email=about['__email__'],
    description=about['__description__'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    url=about['__url__'],
    packages=setuptools.find_packages('pybind11'),
    install_requires=INSTALL_REQUIRES,
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    keywords=about['__keywords__'],
    ext_modules=[relaxation]
)
