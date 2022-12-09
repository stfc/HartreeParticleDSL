#!/usr/bin/env python3
"""Setup script. Used by easy_install and pip."""

import os
from setuptools import setup, find_packages

BASE_PATH = os.path.dirname(os.path.abspath(__file__))
SRC_PATH = os.path.join(BASE_PATH, "src")
PACKAGES = find_packages(where=SRC_PATH)

NAME = 'HartreeParticleDSL'
AUTHOR = ("Aidan Chalk <aidan.chalk@stfc.ac.uk>")
AUTHOR_EMAIL = 'aidan.chalk@stfc.ac.uk'
URL = 'https://github.com/NYI'
DOWNLOAD_URL = 'https://github.com/NYI'
DESCRIPTION = ('HartreeParticleDSL - A Generic Particle DSL supporting a variety of backends')
LONG_DESCRIPTION = '''\
TBD
'''
LICENSE = ' TBD '

CLASSIFIERS = [
    'Development Status :: 3 - Alpha',
    'Environment :: Console',
    'Intended Audience :: Developers',
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Programming Language :: Python',
    'Topic :: Scientific/Engineering',
    'Topic :: Software Development',
    'Topic :: Utilities',
    'Operating System :: POSIX',
    'Operating System :: Unix']

VERSION = '0.0.1a'

if __name__ == '__main__':

    def get_files(directory, install_path, valid_suffixes):
        '''Utility routine that creates a list of 2-tuples, each consisting of
        the target installation directory and a list of files
        (specified relative to the project root directory).

        :param str directory: the directory containing the required files.
        :param str install_path: the location where the files will be placed.
        :param valid_suffixes: the suffixes of the required files.
        :type valid_suffixes: [str]

        :returns: a list of 2-tuples, each consisting of the target \
            installation directory and a list of files (specified relative \
            to the project root directory).
        :rtype: [(str, [str])]

        '''
        examples = []
        for dirpath, _, filenames in os.walk(directory):
            if ("__" not in dirpath) and filenames:
                rel_path = os.path.relpath(dirpath, directory)
                files = []
                for filename in filenames:
                    if any([filename.endswith(suffix) for
                            suffix in valid_suffixes]):
                        files.append(
                            os.path.join(os.path.basename(install_path),
                                         rel_path, filename))
                if files:
                    examples.append((os.path.join(install_path, rel_path),
                                     files))
        return examples

    # We have all of the example, tutorial and wrapper libraries files
    # listed in MANIFEST.in but unless we specify them in the data_files
    # argument of setup() they don't seem to get installed.
    # Since the data_files argument doesn't accept wildcards we have to
    # explicitly list every file we want.
    # INSTALL_PATH controls where the files will be installed.
    # VALID_SUFFIXES controls the type of files to include.

    EGS_DIR = os.path.join(BASE_PATH, "examples")
    INSTALL_PATH = os.path.join("share", "HartreeParticleDSL", "examples")
    VALID_SUFFIXES = ["90", "py", "md", ".c", ".cl", "Makefile", ".mk", ".cpp", ".hpp"]
    EXAMPLES = get_files(EGS_DIR, INSTALL_PATH, VALID_SUFFIXES)

    LIBS_DIR = os.path.join(BASE_PATH, "lib")
    INSTALL_PATH = os.path.join("share", "HartreeParticleDSL", "lib")
    VALID_SUFFIXES = ["90", "sh", "py", "md", "Makefile", ".mk",
                      ".jinja", "doxyfile"]
    LIBS = get_files(LIBS_DIR, INSTALL_PATH, VALID_SUFFIXES)
    setup(
        name=NAME,
        version=VERSION,
        author=AUTHOR,
        author_email=(AUTHOR_EMAIL),
        license=LICENSE,
        url=URL,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        classifiers=CLASSIFIERS,
        packages=PACKAGES,
        package_dir={"": "src"},
        install_requires=['pyparsing', 'six'],
        extras_require={
            'dag': ["graphviz"],
            'doc': ["sphinx", "sphinxcontrib.bibtex < 2.0.0",
                    "sphinx_rtd_theme", "autoapi"],
            'psydata': ["Jinja2"],
            'test': ["pep8", "pylint", "pytest-cov", "pytest-pep8",
                     "pytest-pylint", "pytest-flakes", "pytest-pep257"],
        },
        include_package_data=True,
#        scripts=['bin/psyclone', 'bin/genkernelstub', 'bin/psyad'],
        data_files=LIBS
        )

