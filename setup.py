#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['Click>=7.0', ]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest>=3',]

docs_requirements = [dependency for dependency in open("requirements_docs.txt").readlines()]

dev_requirements = [dependency for dependency in open("requirements_dev.txt").readlines()]

setup(
    author="David Huard",
    author_email='huard.david@ouranos.ca',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="A Python wrapper to setup and run the hydrologic modelling framework Raven.",
    entry_points={
        'console_scripts': [
            'ravenpy=ravenpy.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type="text/x-rst",
    include_package_data=True,
    keywords='ravenpy',
    name='ravenpy',
    packages=find_packages(include=['ravenpy', 'ravenpy.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    extras_require={
        "docs": docs_requirements,
        "dev": dev_requirements,
    },
    url='https://github.com/CSHS-CWRA/ravenpy',
    version='0.1.0',
    zip_safe=False,
)
