import os
import sys

from setuptools import find_packages, setup

import versioneer


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="litebirdtask",
    provides="litebirdtask",
    version=versioneer.get_version(),
    description="LiteBIRD TOAST Analysis and Simulation Kit",
    long_description=readme(),
    long_description_content_type="text/markdown",
    author="LiteBIRD Collaboration",
    author_email="tskisner.public@gmail.com",
    url="https://github.com/hpc4cmb/litebirdtask",
    packages=["litebirdtask"],
    scripts=None,
    license="BSD",
    python_requires=">=3.6.0",
    setup_requires=["wheel"],
    install_requires=["toast-cmb"],
    cmdclass=versioneer.get_cmdclass(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Operating System :: POSIX",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Astronomy",
    ],
)
