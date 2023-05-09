HartreeParticleDSL internals
============================

HartreeParticleDSL has undergone rapid changes during its lifetime, and the current
best supported implementation uses the "Particle IR" intermediate represention.

Particle IR is heavily based on PSyclone's PSyIR.

To construct the intermediate representation, the AST_to_Particle_IR module is used.
This uses Python's `ast` module to walk through kernels, and produce Particle IR nodes
that correspond to that code.

Particle IR is primarily used by HartreeParticleDSL backends to generate code, but
they also do some analysis to do other tasks, e.g. automatically add boundary conditions
into code when relevant or perform MPI communication when needed.

In the future, Particle IR should also enable optimisations, transformations, kernel
merging etc.


Kernel Structure in Particle IR
-------------------------------
Each type of kernel that is supported in Particle IR have a special node in the Particle
IR, which require a specific minumum number of arguments. Backends may also require or allow
extra arguments for some kernels, so check their documentation for a detailed discussion of
the arguments required/allowed for each kernel type.

Each kernel specified in a user's program is converted into the relevant set of Particle IR nodes,
which are then stored in the backend ready for code generation. These kernels are entirely
independent trees, and are usually operated on as separate entities.

In the future, backends may choose to combine/merge or perform other operations on these kernels
to optimise them for the provided implementations, the kernels are left as implemented by the user.

Typing and types in Particle IR
-------------------------------
The Particle IR contains a number of `datatype` implementations, which handles how datatype are
created and stored in the Particle IR. Custom datatype objects are fairly easy to create for
the user, however this not yet documented.

Datatypes are divided into a few base types:

1. NoType (e.g. `void` type in C)
2. ScalarType - This is used for most common scalar datatypes, e.g. `int` or `double`. These datatypes
   have an `Intrinsic` type, e.g. `INTEGER` or `FLOAT` etc. and a `Precision`. Backends then convert
   this representation into the appropriate output for the chosen language.
3. StructureType - Structures are made of a structure type containing a set of components. The
   components are contained in a dict of `(Name: str, DataType)` pairs. StructureTypes can contain
   other structures as components, as is common in most languages.
4. PointerType - This represents a pointer to a datatype. Mostly this is used by backends directly
   if at all, and is primarily used in case compatibility is required.
5. ArrayType - This represents an Array of a DataType, and is made up of a datatype and shape. The
   shape is a list of array extents, with `Extent.DYNAMIC` being available for runtime defined extents.
   Array indices are internally expected to start from 0 and run to extent-1 for each dimension, using
   a row-major ordering (i.e. C-style).

As well as the base types, Particle IR defines a set of inbuilt scalar datatypes which should be sufficient
for most use-cases.

The Particle and Config types in the main DSL have their own types, which have helper functions
to add elements to. For common C/C++ types these can then be added to the particle or config
easily without needing to involve Particle IR.

Symbols in Particle IR
----------------------
Symbols are divided into types, which mostly correspond to the available datatypes, e.g. a variable
with a `ScalarType` will be a `ScalarTypeSymbol`. There is one additional symbol type without a
corresponding datatype, which is `AutoSymbol`. This is used to provide an interface with complex
C++ calls, usually which are built into a backend as opposed to being constructured from user code.
One example of its use case is to support random number generation in kernels inside the Cabana 
backend, as supporting this needs `auto generator = random_pool.get_state()`, as the resulting type
is too complex otherwise to easily represent.

AST_to_Particle_IR module
-------------------------
This module takes each module from the Python AST layout, and converts it into Particle IR.
Not all of Python's functionality is supported by HartreeParticleDSL, and unsupported statements
will throw an exception from the `generic_visit` function.

Most supported nodes have an obvious transformation to Particle IR, however this document notes a few
which have more complicated transformations.

If Statements
^^^^^^^^^^^^^
In Python, if statements can be written with multiple `elif` conditions followed by an `else` condition.
In Particle IR, these statements will be stored as multiple nested `IfElseBlock` nodes, where the
first statement in the else block is the if condition from the `elif` test.

Call Statements
^^^^^^^^^^^^^^^
HartreeParticleDSL has a few special call types that are reserved, and the Particle IR construction handles
these separately.

The first is the `create_variable` call. In this case, the variable is added to the containing kernel/function's
symbol table, and an `EmptyStatement` node is added to the tree. In most cases, this statement can then be
ignored during code generation. If the `create_variable` call has an initial value assignment, the node is
instead replaced by an `Assignment` node to set the value of the created variable at this point in the tree.

Next is the `invoke` call. This maps straight to the `Invoke` member in the Particle IR, and is then handled
during codegen.

The final special call is `initialise`. This call takes key word arguments so needs to be handled differently,
however it does produce a standard `Call` node.

Backend structure
-----------------

Backends that support Particle IR need a PIR_Visitor subclass, which handles the conversion of
Particle IR trees into the relevant code. This will usually contain visitors for almost all
nodes supported by Particle IR, though it doesn't have to.

Cabana_PIR_backend
^^^^^^^^^^^^^^^^^^
This is currently the only non-deprecated backend, as it is the only backend to support PIR.
The backend supports all PIR nodes except `PairwiseKernel`.
The output from this backend is C++ code compatible with Cabana (https://github.com/ECP-copa/Cabana)
and Kokkos (https://github.com/kokkos/kokkos). It can output code that supports both CUDA GPUs and
MPI parallelisation.

It has several inbuilt functions such as `get_random_number`, `pi`, `println`.

It supports Parallel HDF5 for I/O.
