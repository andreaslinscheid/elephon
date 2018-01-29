# Elephon

Computing electron-phonon interactions from first principles.
The code uses finite displacements to compute all phonon-related properties, including the linear
displacement potential created by a phonon mode. For details on the underlying method, see [to be published].

The goal of this software is to compute the electron-phonon superconducting critical temperature.
Initially, this will be done by computing the isotropic electron-phonon Eliashberg function.

PLEASE NOTE: This is a prelinmary development version.
In this initial stage, we establish that (and how) the principle method works and the focus is not on maintainability.
In the next stage, we will extend the usability of the software and refactor to achieve high maintainability.

## Project goals
* Provide a reliable and numerically efficient postprocessing software to compute electron-phonon superconducting critical temperature
* Work with [VASP](https://www.vasp.at/) as the electron-structure calculator.
* Based on an initial ground state calculation, set-up input files to make the calculation with minimal user interaction.

## Design principles
* Since the extent of the calculation is significant, we require performant compiled code (C++).
  The focus will initially be on a linux operating system using x86 architecture.
* Easy to build and install for the user (build process using [CMake](https://cmake.org/))
* Maximal defaults: Try to specify default input parameters so the user has to configure as little as possible while still keeping the software flexible.
* Provide easy to use and rich graphical ouput methods (plots of displacment potentials ect) to allow experienced users to understand the physics of a system.
* Code documentation in-place from annotated source code via [Doxygen](www.doxygen.org/).
* Revision controll using git.

## Prerequisites
We link to the [BOOST](www.boost.org/) software libraries, version 1.56 or greater. Furthermore, 
this package requires [FFTW3](www.fftw.org/) and BLAS/LAPACK implementations.
For graphical output, we require the [VTK](https://www.vtk.org/) libraries.
