# Automated testing and making template examples (beta version)
**One can unit test mechanoChemFEM and obtain some working examples in 2 easy steps**

**Step 1**: Pulling unit test scripts and data from GitHub Repo
 - Using svn (install in Mac: sudo port install subversion) to only download unit-tests subfolder from GitHub
```
svn checkout https://github.com/mechanoChem/morphogenic-interfaces/trunk/unit-tests
```
 - Clone entire repo including unit test subfolderCancel changes
```
git clone https://github.com/mechanoChem/morphogenic-interfaces.git
```

**Step 2**: runnig the `unit_test.py` command inside the downloaded unit-tests folder
```
python unit_test.py my_test_folder test-1
```
- unit_test.py has 2 inputs (more options as inputs will be added)
- 1st  input: path/to/test/folder (here test_folder), where the unit test files and results are generated 
- 2nd input: test name, currently test-1 is available

# Dependencies
[dealii](https://www.dealii.org/), [mechanoChemFEM](../lib), python 3.8+ and numpy, [meshio](https://github.com/nschloe/meshio) (pip install meshio, MIT license) 




