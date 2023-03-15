# Extended-OP2
Extending OP2 Framework to Support Portable Parallel Programming of Complex Applications

The original OP2 and more details can see in https://github.com/OP-DSL/OP2-Common

## Modify
We extend the source to source translator and the runtime library of OP2,which makes it possible to automatically support these complex applications with complex data structures such as multidimensional arrays and structures and complex function calls in the generation of serial, OpenMP, CUDA and MPI versions of codes.

## Start
Firstly, refer to https://github.com/OP-DSL/OP2-Common to install OP2.

Then, Using OP2 to automatically generate parallel code(OpenMP, MPI and CUDA) of complex application.

Finally, to build OP2 again and application.

## Citing
To cite OP2, please reference the following paper:

G. R. Mudalige, M. B. Giles, I. Reguly, C. Bertolli and P. H. J. Kelly, "OP2: An active library framework for solving unstructured mesh-based applications on multi-core and many-core architectures," 2012 Innovative Parallel Computing (InPar), 2012, pp. 1-12, doi: 10.1109/InPar.2012.6339594.
