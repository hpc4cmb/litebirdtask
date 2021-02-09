# LiteBIRD TOAST Analysis and Simulation Kit

This package provides tools for simulation and analysis of LiteBIRD data within the
[TOAST](https://github.com/hpc4cmb/toast) framework.

## Required Inputs

This repository contains no private information about LiteBIRD.  Before using these
tools you must [download a model of the instrument from this page](https://wiki.kek.jp/display/cmb)

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

Next, use pip to install this package and its requirements:

    pip install litebirdtask

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
