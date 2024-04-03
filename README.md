# OrthoBase

A Python library for constructing multiplet projectors in the decomposition of
the product of adjoint representations in the SU(N_c) group. The library also
provides numerous methods to work with Young tableaux. The software achieves
high efficiency through a combination of the
[Z3](https://www.microsoft.com/en-us/research/project/z3-3/) satisfiability
modulo theories solver and the symbolic manipulation program
[FORM](https://www.nikhef.nl/~form/).

Key features:
- Construct multiplet projectors for the decomposition of adjoint
representation products in SU(N_c) groups
- Efficient algorithms implemented using the Z3 solver and FORM symbolic
manipulation
- Comprehensive set of methods for working with Young tableaux
- Easy-to-use Python interface for seamless integration into your projects
- Well-documented code and extensive online documentation for quick start and
reference

The library's toolset provides researchers working in fields involving SU(N_c)
groups and Young tableaux with a powerful and efficient means to streamline
their work in these domains.

# Installation

OrthoBase requires Python version 3.10 or newer and can be installed on any
Unix-like operating system.

The latest release can be installed from [PyPI] using [pip]:

	$ python3 -m pip install --user --upgrade OrthoBase

This will install the library with all its python dependencies. There are
additional dependencies described below. They should be installed manually or
using the package manager of your operating system.

[pip]: https://pypi.org/project/pip/
[pypi]: https://pypi.org/project/OrthoBase/

## Additional dependencies

These dependencies are not automatically installed via pip. Therefore, users
must manually install them and ensure that the relevant files are included in
the `$PATH` and `$LD_LIBRARY_PATH` environment variables:

- FORM (http://www.nikhef.nl/~form/)
- An implementation of MPI (Message Passing Interface) should be installed. Some popular implementations include:
  - [Open MPI](https://www.open-mpi.org/)
  - [MPICH](https://www.mpich.org/)
  - [Intel MPI](https://software.intel.com/en-us/mpi-library)
  - [Microsoft MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi)

# Usage

Refer to the [online documentation](orthobase.readthedocs.io) of our package.

You can run the example calculation as follows:

	$ ./demo.py
