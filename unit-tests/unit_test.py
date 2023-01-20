#!/usr/bin/env python3

'''
Unit testing mechanochemFEM (beta version)
Mostafa Shojaei, Jan 2023

Dependencies:
python 3.8+, numpy, meshio (pip install meshio)
'''

# importing python modules
import os
import sys
import shutil
import subprocess
import numpy as np
import re
import time
import meshio


def run_shell(shell_command,output_file_path=None):
    if output_file_path == None:
        subprocess.call(shell_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True)
    else:
        output_file = output_file_path.rsplit("/",1)[1]
        with open(output_file_path, "a") as file:
            subprocess.call(shell_command + " >> " + output_file,
                stdout=file,
                stderr=file,
                shell=True)
    
def run_shell_1(shell_command,output_file_path=None):
    if output_file_path == None:
        return subprocess.Popen(shell_command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        shell=True)
    else:
        output_file = output_file_path.rsplit("/",1)[1]
        with open(output_file_path, "a") as file:
            return subprocess.Popen(shell_command + " >> " + output_file,
                stdout=file,
                stderr=file,
                shell=True)

def print_process(process,msg,delay=0.5):
    i = '   '
    while process.poll() is None:
        print(f' {msg}{i}\r', end='')
        if    i == '   ':  i = '.  '
        elif  i == '.  ':  i = '.. '
        elif  i == '.. ':  i = '...'
        else: i ='   '
        time.sleep(delay)
    print(f' {msg}(done)')

try:
    dir = sys.argv[1]
except:
    print('No input dir.\n')
else:
    print('\nCode directory:',dir)
    try:
        test_0 = sys.argv[2]
    except:
        test_0 = 'test_0'

    ref_folder = 'ref'

    copy = True
    comp = True
    run = True
    show_run_prgress = True
    verify = True

    if len(sys.argv[0].rsplit("/",1)) == 2:
        ref_folder = sys.argv[0].rsplit("/",1)[0] + "/" + ref_folder

    vtk_folder = f'{ref_folder}/{test_0}/vtk_results'
    build_dir = dir + '/build'
    # compile_outputs = build_dir + '/compile_outputs.txt'

    # making a folder to run the tests
    print(f'\nSetting up {dir}')
    if not os.path.exists(dir):
        print(f"  - {dir} made")
        os.mkdir(dir)
    else:
        print(f"  - {dir} already exists")
    if not os.path.exists(build_dir):
        print(f"  - {build_dir} made")
        os.mkdir(build_dir)
    else:
        print(f"  - {build_dir} already exists")

    # copying required files from ref folder
    if copy:
        print(f'  - code scripts copied from the reference folder to {dir}',end='')
        shutil.copyfile(f'{ref_folder}/{test_0}/main.cc',f'{dir}/main.cc')
        shutil.copyfile(f'{ref_folder}/{test_0}/CahnHilliard.h',f'{dir}/CahnHilliard.h')
        shutil.copyfile(f'{ref_folder}/{test_0}/CMakeLists.txt',f'{build_dir}/CMakeLists.txt')
        shutil.copyfile(f'{ref_folder}/{test_0}/parameters.prm',f'{build_dir}/parameters.prm')
        print(' (done)')
        # with open(f'{build_dir}/parameters.prm','r') as pf:
        #     print(pf.read())

    # compile the code
    if comp:
        print(f'\nCompiling the code in {build_dir}')
        p = run_shell_1(f'cd {build_dir}; cmake CMakeLists.txt', f"{build_dir}/cmake_compile_outputs.txt")
        print_process(p," - Generating make files using cmake ")
        p = run_shell_1(f'cd {build_dir}; make release; make;', f"{build_dir}/make_compile_outputs.txt")
        print_process(p," - Running make to compile the code ")
        
    # run
    if run:
        print(f'\nRunning the code in {build_dir} ...')
        test_output_path = f"{build_dir}/0_outputs.txt"
        print(f'  - See outputs in {test_output_path} (refresh to see updates)')
        lines = 20
        print(f'  - Use following in a seprate terminal to see last {lines} lines of output:')
        print(f'    watch tail -n {lines} {os.getcwd()}/{test_output_path}')
        p = run_shell_1(f'cd {build_dir}; mpirun -np 32 main -pc_type lu -pc_factor_mat_solver_type superlu_dist', test_output_path)
        if show_run_prgress:
            vtk_files = os.listdir(vtk_folder)
            vtk_files.sort(key=lambda f: int(float(re.sub('\D', '', f))))
            n = int(float(re.sub('\D', '', vtk_files[-1])))
            if not os.path.exists(build_dir+"/vtk_outputs/"):
                os.mkdir(build_dir+"/vtk_outputs/")
            k = 0
            while k < n:
                print(f'  - Progress {k/n*100:.2f}%\r', end='')
                # sys.stdout.write('\r')
                # sys.stdout.write("%d%%" % (k/n*100))
                # sys.stdout.flush()
                k = len(os.listdir(build_dir+"/vtk_outputs/"))
            print(f'  - Progress {100.00:.2f}%\r', end='')
        else:
             print_process(p," - this may take several hours ")

    # verify
    er = 1e-6
    rmse = lambda x: np.sqrt(np.sum(x**2))
    passed = True
    if verify:
        print(f'\nVerifying .vtk results generated in {vtk_folder} ...')
        vtk_files = os.listdir(vtk_folder)
        vtk_files.sort(key=lambda f: int(float(re.sub('\D', '', f))))
        for vtk_file in vtk_files:
            
            ref_path = vtk_folder+"/"+vtk_file
            res_path = build_dir+"/vtk_outputs/"+vtk_file
            print(f'Checking vtk result {res_path}',end=', ')
            if os.path.exists(res_path):
                # load ref
                vtk_data_ref = meshio.read(ref_path)
                points_ref = vtk_data_ref.points
                values_ref = vtk_data_ref.point_data
                x_ref = points_ref[:,0]
                y_ref = points_ref[:,1]
                c_ref = values_ref['c1'].T[0]
                # pc_ref = vtk_data_ref.cells_dict['quad'].T  # ponints_elements connectivity table
                
                # load result
                try:
                    vtk_data = meshio.read(res_path)
                except:
                    passed = False
                    print('(failed)')
                    print(f'Generated result {res_path} is not readable')
                    break
                else:
                    points = vtk_data.points
                    values = vtk_data.point_data
                    x = points[:,0]
                    y = points[:,1]
                    c = values['c1'].T[0]
                    # pc = vtk_data.cells_dict['quad'].T  # ponints_elements connectivity 
                    
                    if rmse(x-x_ref)==0 and rmse(y-y_ref)==0:
                        if rmse(c-c_ref) < er:
                            pass
                        else:
                            passed = False
                            print('(failed)')
                            print(rmse(c-c_ref))
                            print(f'{res_path} is different from its reference {ref_path}')
                            break
                    else:
                        passed = False
                        print('(failed)')
                        print(f'mesh in {res_path} is different from its reference {ref_path}')
                        break
                print('(passed)')
            else:
                passed = False
                print('(failed)')
                print(f'Expected result {res_path} does not exist')
                break
        
        if passed:
            print('------\nTest passed\n------')
        else:
            print('------\nTest failed\n------')
            



                
