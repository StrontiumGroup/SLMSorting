.. pySLMSorting documentation master file, created by
   sphinx-quickstart on Fri Apr 18 12:11:58 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pySLMSorting documentation
==========================

This is the documentation for the python version of the SLMSorting code made at the University of Amsterdam. See also the publication at `ArXiv <https://arxiv.org/abs/2501.01391>`_.

While the main code was written in C++, a simple debugging and testing code was written in python. Note that this does not claim nearly the same performance as the C++ code!

This documentation covers three parts:

* An overview of the functions used in the python module;
* A description of the kernels that are used in the GPU;
* An example notebook.

Note that unlike the C++ code, the python code does not include some class that automatically does the linear interpolation. This is because at the UvA, we used the python code mostly for calculating single holograms using the WGS algorithm. A simple example on how to use the kernels to do the interpolation is presented in the code.

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   :glob:
   
   modules
   kernels
   Example


