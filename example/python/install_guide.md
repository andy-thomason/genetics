# Boost.Genetics Python Wrapper
The enclosed bindings.cpp file will use Boost Python to wrap Boost.Genetics C++
 classes into Python classes and bundle them into a single importable Python module.

This document aims to describe how to build and integrate this Boost.Genetics into
a larger Python package.

## Installation

Installation is divided into two parts. First, our goal is to convert bindings.cpp
 into a dynamically-loadable Python module called "genetics" 
 (.dll, .so, .dylib depending on your OS).
While this is not a "typical" .py module, basic usage should be as simple as:

```python
cd example/python
make
python
import genetics
```
Unlike Python C extensions, Cython, CFFI, or other Python-driven, Boost.Python
is a C++ driven code **generation** tool. Thus, this first goal will largely be
driven by the C++ tool chain and dependencies. The largest of these is a C++
compiler and Boost itself.

### Dependencies
- Boost.Python
- Boost.Genetics
- Python
- The Clang C++ compiler 


### Python integration and distribution

Some users may want to go further than simply using Boost.Genetics inside of a
Python script by incorporating it into a larger Python package. To do this, you
will need to switch from the C++ tool chain to the Python tool chain.

In this phase, our goal is to take the Boost-generated genetics.so module and 
install it into our PYTHONPATH so it may be imported by other parts
of a Python project.

This may be done simply by moving the genetics.so file around, but such a
solution will make it difficult to remove, upgrade, re-install, and otherwise 
manage future releases.

A second, simple option is to distribute pre-built genetics.so/dll/pyd/dylib 
binaries for the target platform.

For the brave (or developers), a third option is to write a setup.py script with
a custom install command. This setup script can live in a shell python package
and then be invoked using familiar PIP syntax (pip install boost-genetics-py)

The general idea is to allow a user to specify a python package as a requirement
and automate the build and 

1) Clone a release of the Boost.Genetics C++ code and mark the location of the
"include/boost/genetics" directory.

2) Assert a suitable build environment is present. This is very platform dependent.

3) Generate a "Pythonic" directory structure by linking bindings.cpp and test.py
to their appropriate locations.

|- setup.py
|- extensions
|  |- bindings.cpp
|- test.py
  
4) Configure the package meta data (version of Boost.Genetics) and any other
build parameters.

5) Compile the extension with setup.py build

6) Allow setuptools to do the work of registering the genetics module on the 
appropriate PYTHONPATH

7) Run a basic python test to make sure everything is sane 

##Building with b2

This project should live in your_boost/libs/genetics

The Jamfile.v2 file should be enough to build the DLL

On windows:

..\..\..\..\..\b2 address-model=64

This should make genetics.pyd

*note* it is essential that

a) You use a 64 bit Python build
b) You run the bootstrap first and build the boost libs
c) You make sure the boost python DLL is in the path

If python complains about missing DLLs, it is probably not your
DLL but the boost DLLs which are inacessible.

Once it works, test with "import genetics" in the python
command shell.



