<B>mechanoChemFEM</B><br>
=======================================================================
mechanoChemFEM is a libarary for modeling of mechano-chemical problems using the finite element method. It consists of [classes](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/annotated.html) and [functions](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/modules.html) based on [Deal.ii](https://www.dealii.org/), and examples.

<B>Master branch</B> contains the current stable version of the code and documentaion.<br>
<B>example branch</B> contains all example code.<br>
<B>Other branches</B>  contain the older versions.


<B>List of contributors:</B><br>

Zhenlin Wang (lead developer）<br>

Krishna Garikipati<br>

[Computational Physics Group, University of Michigan](http://umich.edu/~compphys/index.html)


[<B>Code documentation</B>](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/index.html)

<B>Overview</B><br>
=======================================================================
mechanoChemFEM lib contains six modules:


	hpFEM : offer high level functions dedicated to mulitphysics modeling using deal.ii
	
	Residual : offer residual functions of a variety of equations and boundary conditions

	solveClass : offer nonlinear solvers based on deal.ii linear solver

	FEMdata : high level functions for data output and restart(checkpoint)
	
	InitBoundValProbs :  a pre-defined generic initial boundary value problem framework
	
	supplementary :a collection of extra data structures and many helper functions

	
<B>ParameterHandler</B> is used for parameter management. 


#  Installation
1. Install dependencies
	  1) Install [CMake](http://www.cmake.org/download/)
	  2) Install [deal.II](www.dealii.org/download.html) with [Trilinos](https://trilinos.org/) and [PetSc](https://www.mcs.anl.gov/petsc/download/index.html) (Deal.II OSX binaries include full packages of deal.ii with Trillions and other useful libs.)
2. Install mechanoChemFEM
	  1) cd into “build” folder
	  2) Modify CMakeList.txt for path of pre-required libs: deal.ii (with Trilinos, Petsc)
	  3) cmake CMakeLists.txt
	  4) `make install` or do `make release install`
	  5) `make run`
	  
# Usage
Following examples discuss in details how to implement PDEs, set input parameters, and obtain the results in mechanoChemFEM\
[Example 1: Diffusion-reaction equation](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/diffusion_reaction.html)\
[Example 2: Cahn-Hilliard equation](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/_cahn_hilliard.html)\
[Example 3: Allen-Cahn equation](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/_allen__cahn.html)\
[Example 4: Coupled system with multiple domains](https://htmlpreview.github.io/?https://raw.githubusercontent.com/mechanoChem/mechanoChemFEM/master/doxygen/html/growth.html)\

<B>License</B><br>
=======================================================================
GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.

<B>Acknowledgements</B><br>
=======================================================================
This code has been developed under the support of the following: <br>

- NSF DMREF grant: DMR1436154 "DMREF: Integrated Computational Framework for Designing Dynamically Controlled Alloy-Oxide Heterostructures" <br>
- DOE BES, Division of Materials Sciences and Engineering: Award #DE-SC0008637 that funds the PRedictive Integrated Structural Materials Science (PRISMS) Center at University of Michigan <br>
- Toyota Research Institute, Award #849910, "Computational framework for data-driven, predictive, multi-scale and multi-physics modeling of battery materials" <br>
