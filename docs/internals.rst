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

The Particle and Config types in the main DSL have their own types, which have helper functions
to add elements to. For common C/C++ types these can then be added to the particle or config
easily without needing to involve Particle IR.

Backend structure
-----------------

Backends that support Particle IR need a PIR_Visitor subclass, which handles the conversion of
Particle IR trees into the relevant code. This will usually contain visitors for almost all
nodes supported by Particle IR, though it doesn't have to.
