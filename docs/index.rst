.. HartreeParticleDSL documentation master file, created by
   sphinx-quickstart on Sat Oct  2 04:03:45 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to HartreeParticleDSL's documentation!
==============================================

HartreeParticleDSL is a particle DSL written in Python. The goal of the library
is to provide a simple high level front-end for writing particle codes, whilst
enabling access to high-performance parallel computers. This is done by providing
backends to target high performance runtime libraries, e.g. FDPS. 

The basic premise behind describing particle codes comes through describing the 
kernels used in the computation. Currently HartreeParticleDSL supports:

1. Pairwise Interactions. These are interactions that occur over particle pairs 
   within a cutoff radius.
2. Per-particle kernels. These are kernels executed on every particle in the
   system. A common-use case for this type of kernel is for timestepping loops.
  

Examples are found in the examples folder. These are currently relatively disorganised,
and in future we will expand these to better document simple use-cases for each backend
and IO module.


Each backend can have their own support of python language features, however as
much code as possible should be portable between backends. One example of how this
portability may not exist could be between a simple C and Fortran backend. A line
of Python code:

.. code-block:: Python

  array[:] = array[:] * 2.0

Since Fortran has language support for array notation, a Fortran backend could choose
to support Python's array-slice functionality. For a C-based backend to support such
notation would require significantly more information, so may choose not to.

More generic language features (such as for/while loops, if statements, multidimensional
arrays) must be supported by all backends - the requirements for backend support is
currently undocumented, but will be added in future.




.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
