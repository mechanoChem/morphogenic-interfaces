# morphogenic-interfaces
This repository provides **phase-field simulations** based on [mechanoChemFEM](https://github.com/mechanoChem/mechanoChemFEM)
for design and study of intermetallic solid/solid charge transfer interface materials to improve the performance of Li-ion solid-state batteries

Main directories: 
- [unit-tests](unit-tests/): Automated framework to test mechanoChemFEM and provide template codes for further development by user
- [case-studies](case-studies/): Codes and input parameters for all case studies 
- [lib](lib/): [mechanoChemFEM](https://github.com/mechanoChem/mechanoChemFEM) latest source code compatible with all studies/tests
- [utilty](utilty/): useful scripts for pre/post-processing 

## mechanoChemFEM
mechanoChemFEM is a comprehensive open source library (licensed by LGPL) for modeling of mechano-chemical problems using the finite element method developed in [Computational Physics Group at University of Michigan](http://umich.edu/~compphys/index.html), and It consists of [classes](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/annotated.html) and [functions](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/modules.html) based on [Deal.ii](https://www.dealii.org/). 

- **A detailed documentation of mechanoChemFEM is available** [here](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/index.html)>
 - **The latest version of mechanoChemFEM compatible with problems studied in this repo is included in [lib](lib/) directory**. 

###  Installation
1. Install dependencies
	  1) Install [CMake](http://www.cmake.org/download/)
	  2) Install [deal.II](www.dealii.org/download.html) with [Trilinos](https://trilinos.org/) and [PetSc](https://www.mcs.anl.gov/petsc/download/index.html) (Deal.II OSX binaries include full packages of deal.ii with Trillions and other useful libs.)
2. Install mechanoChemFEM
	  1) cd into “build” folder
	  2) Modify CMakeList.txt for path of pre-required libs: deal.ii (with Trilinos, Petsc)
	  3) cmake CMakeLists.txt
	  4) `make install` or do `make release install`
	  5) `make run`

### Usage
Following examples discuss in details how to implement PDEs, set input parameters, and obtain the results in mechanoChemFEM\
[Example 1: Diffusion-reaction equation](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/diffusion_reaction.html)\
[Example 2: Cahn-Hilliard equation](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/_cahn_hilliard.html)\
[Example 3: Allen-Cahn equation](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/_allen__cahn.html)\
[Example 4: Coupled system with multiple domains](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/growth.html)\
