# LiteBIRD TOAST Analysis and Simulation Kit

**WARNING:  THIS PACKAGE IS NOT YET RELEASED OR READY.**  There will be communication to the mailing lists and [updates on the wiki page](https://wiki.kek.jp/pages/viewpage.action?pageId=150667506) when it is.

This package provides tools for simulation and analysis of LiteBIRD data within the
[TOAST](https://github.com/hpc4cmb/toast) framework.

## Required Inputs

This repository contains no private information about LiteBIRD.  You should [visit the
official wiki page here](https://wiki.kek.jp/pages/viewpage.action?pageId=150667506) and download the required
files for use with this package.

## Installation

First, make sure that you have a recent (>=3.6) version of Python3.  This might come
from your operating system, a package manager, or Anaconda.  You can check your python version with:

    python3 --version

Next create a virtualenv (name it whatever you like):

    python3 -m venv ${HOME}/litebird

Now activate this environment:

    source ${HOME}/litebird/bin/activate

Within this virtualenv, update pip to the latest version:

    python3 -m pip install --upgrade pip

Now you can choose whether to install the latest stable version of the package from
PyPI, or to install from a git checkout.  To install a stable version, do:

    pip install litebirdtask

*OR*, from your git checkout of `litebirdtask`, do:

    pip install .

If you wish to use other tools in this environment, you can install them now with pip.
To "unload" this environment just do:

    deactivate

The next time you wish to use the tools, load them again:

    source ${HOME}/litebird/bin/activate

You can always completely delete this environment by simply deleting the directory you
created.  There is extensive documentation on the internet about managing python virtual
environments.

## Documentation

The full documentation is [available on readthedocs]().
