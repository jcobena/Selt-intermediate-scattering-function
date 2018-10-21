'''
Created on Nov 1, 2017

@author: Jose_Cobena
'''
import sys, csv
import numpy as np
import pandas as pd
import itertools
import math as mt
import collections
from itertools import count, zip_longest
#from pandas.core.frame import DataFrame

def main():
    T = 343
    
    initial_step_2 = 0 #1230000
    initial_step_10 = 20
    initial_step_100 = 200
    initial_step_1000 = 2000
    initial_step_10000 = 10000
    initial_step_100000 = 200000
    initial_step_1000000 = 2000000
    initial_step_10000000 = 7000000
    
    final_step_2 = 11
    final_step_10 = 101
    final_step_100 = 1001
    final_step_1000 = 10001
    final_step_10000 = 100001
    final_step_100000 = 1000001
    final_step_1000000 = 6000001
    final_step_10000000 = 36000001
                          

    com_file_2 = 'com{}-2.txt'.format(T)
    com_file_10 = 'com{}-10.txt'.format(T)
    com_file_100 = 'com{}-100.txt'.format(T)
    com_file_1000 = 'com{}-1000.txt'.format(T)
    com_file_10000 = 'com{}-a-2000.txt'.format(T) #valid until 6million
    com_file_7mill = 'com{}-1million.txt'.format(T)
    
    increments_2 = 2
    increments_10 = 10
    increments_100 = 100
    increments_1000 = 1000
    increments_10000 = 10000
    increments_100000 = 100000
    increments_1000000 = 1000000
    #increments_10000000 = 1000000
    
    #cut_step = 10000
    #final_step = 22000#5229000 
    #increments = 2000
    n_atoms = 1702 #hydrogen atoms
    
    Lx = 11.91
    Ly = 11.91
    Lz = 334.8
    L = Lz
 
    choices_of_wn = 4
    
#     xx = np.logspace(1, 10, 50, False, 10)
#     yy = np.geomspace(1, 1000, num = 15, endpoint = True)
#     print(*yy)
    
#########################################################################################################
                
    with open('hyd-msd-and-isf-{:08}.txt'.format(initial_step_2), 'r') as file_t0, \
         open('isf_self_fixed_time_all.txt','w') as isf_results_fixed_time, \
         open('isf_self_fixed_wave_all.txt','w') as isf_results_fixed_wave, \
         open(com_file_2, 'r') as com_file2, \
         open(com_file_10, 'r') as com_file10, \
         open(com_file_100, 'r') as com_file100, \
         open(com_file_1000, 'r') as com_file1000, \
         open(com_file_10000, 'r') as com_file2000, \
         open(com_file_7mill, 'r') as com_file7mill:
        
        #print(com_file2000)
        
        com_list2 = positions_list_com(com_file2)
        com_list10 = positions_list_com7(com_file10)
        com_list100 = positions_list_com7(com_file100)
        com_list1000 = positions_list_com7(com_file1000)
        com_list2000_by_lammps = positions_list_com(com_file2000)
        com_list7mill = positions_list_com0(com_file7mill)
        print(com_list2000_by_lammps)
        com_list10000 = positions_list_com_deleting_to_get_log_values(com_list2000_by_lammps, 10000, 100001, 10000)
        com_list100000 = positions_list_com_deleting_to_get_log_values(com_list2000_by_lammps, 200000, 1000001, 100000)
        com_list1000000s = positions_list_com_deleting_to_get_log_values(com_list2000_by_lammps, 2000000, 6000001, 1000000)
        com_list1000000 = positions_list_com_deleting_to_get_log_values(com_list7mill, 7000000, 36000001, 1000000)
        #com_list10000000 = positions_list_com_deleting_to_get_log_values(com_list2000_by_lammps, 20000000, 36000001, 10000000)
        print(com_list2)
        print(com_list10)
        print(com_list100)
        print(com_list1000)
        print(com_list7mill)
        print(com_list10000)
        print(com_list100000)
        print(com_list1000000s)
        print(com_list1000000)
        com_t0 = com_list2[0][1]
        r0_list = positions_list(file_t0)
        
        data_2 = isf_1D_inside_nanotube(initial_step_2, final_step_2, increments_2, com_list2,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        data_10 = isf_1D_inside_nanotube(initial_step_10, final_step_10, increments_10, com_list10,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        data_100 = isf_1D_inside_nanotube(initial_step_100, final_step_100, increments_100, com_list100,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        data_1000 = isf_1D_inside_nanotube(initial_step_1000, final_step_1000, increments_1000, com_list1000,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        data_10k = isf_1D_inside_nanotube(initial_step_10000, final_step_10000, increments_10000, com_list10000,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        data_100k = isf_1D_inside_nanotube(initial_step_100000, final_step_100000, increments_100000, com_list100000,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        data_6m   = isf_1D_inside_nanotube(initial_step_1000000, final_step_1000000, increments_1000000, com_list1000000s,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        data_1000k = isf_1D_inside_nanotube(initial_step_10000000, final_step_10000000, increments_1000000, com_list1000000,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        #data_10M = isf_1D_inside_nanotube(initial_step_10000000, final_step_10000000, increments_10000000, com_list10000000,com_t0, r0_list, 4 ,Lx, Ly, Lz, n_atoms, 'time')
        
       
        #pd.set_option('display.max_columns', None)
        #pd.set_printoptions(max_rows=20000, max_columns=10000)
        
        pd.set_option('max_columns',None)
        pd.set_option('max_rows',None)
        pd.set_option('display.expand_frame_repr', False)
        
        #for x in data_2:
            #y = x.split()
        #print(x, file = isf_results_fixed_time)
        
        print(data_2, file = isf_results_fixed_time)
        print(data_10, file = isf_results_fixed_time)
        print(data_100, file = isf_results_fixed_time)
        print(data_1000, file = isf_results_fixed_time)
        print(data_10k, file = isf_results_fixed_time)
        print(data_100k, file = isf_results_fixed_time)
        print(data_6m, file = isf_results_fixed_time)
        print(data_1000k, file = isf_results_fixed_time)
        #print(data_10M, file = isf_results_fixed_time)
#-----------------------------------------------------------------------------------------------------------------------------#
           

def isf_1D_inside_nanotube(initial_step, final_step, increments, com_list_x, com_t0, r0_list, choices_of_wn,Lx, Ly, Lz, n_atoms, option_output):
    #l, h, m are the wave numbers
    h = 0
    m = 0
    #make 0 if you want to ignore a dimension
    i = 0 #x
    j = 0 #y
    p = 1 #z
    L = Lz
    kc = 2*np.pi
    #dic_isf_it_k1 = {}
    kc_l = 2*np.pi/L
    
    if choices_of_wn == 1:
        wave_numbers = np.arange(mt.pi, Lz/2.0, mt.pi)
    elif choices_of_wn == 2:
        wave_numbers = np.arange(20.45, 40.91, 20.45)
    elif choices_of_wn == 3:
        wave_numbers = np.arange(20.53, 123.5, 20.53) #based in RD
    elif choices_of_wn == 4:
        wave_numbers = np.array([23, 34, 45, 57, 68, 79, 93, 108])
    elif choices_of_wn == 5:
        wave_numbers = np.arange(0, Lz/2.0, 1)
    
    dic_isf_it_k_1000 = {}
    for step, time_step_com in zip(range(initial_step, final_step, increments),com_list_x):
        filesin = 'hyd-msd-and-isf-%.8d.txt' % step
        #filesin = 'msd-%.8d.txt' % step
        files_tt = open(filesin, 'r')#.readlines()
        keys10 = int(step)
        rt_list = positions_list(files_tt)
        results_t_k = []
        
        #substract the displacement of the center of mass
        com_tt = time_step_com[1]
        #print(com_tt)
        com_diplacement_x = com_tt[0] - com_t0[0]
        com_diplacement_y = com_tt[1] - com_t0[1]
        com_diplacement_z = com_tt[2] - com_t0[2]
       
        if time_step_com[0] == step:
            #print(time_step_com[0])
            for l in wave_numbers:
                isf_i_k = 0.0
                for rj_t0, rj_tt, in zip_longest(r0_list, rt_list):
                    if rj_t0[0] == rj_tt[0]:
                        delta_r = []
                        delta_x = rj_tt[1][0] - rj_t0[1][0] - com_diplacement_x
                        delta_y = rj_tt[1][1] - rj_t0[1][1] - com_diplacement_y
                        delta_z = rj_tt[1][2] - rj_t0[1][2] - com_diplacement_z
                        delta_r.extend((delta_x, delta_y, delta_z))
                        
                        isf_i_k += mt.cos(kc * (h * i * (delta_x)/Lx  +  m * j * (delta_y)/Ly  + l * p * (delta_z)/Lz))
                    #print(isf_i_k)
                isf_ave_k = isf_i_k/n_atoms
                results_t_k.append(isf_ave_k)
            #(results_t_k)
            dic_isf_it_k_1000[keys10] = results_t_k
        else:
            #break
            print('there is an error')
    #print(dic_isf_it_k_2)
    od = collections.OrderedDict(sorted(dic_isf_it_k_1000.items()))
    data_from_dict_time_in_x = pd.DataFrame.from_dict(od, orient = 'index')
    #data_from_dict_time_in_x.index = relevant_peaks
    #data_from_dict_time_in_x.index.name = 'time-step'
      
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    
    
    if option_output == None or option_output == 'time':
        #print(data_from_dict_time_in_x.to_string(), file = isf_results_fixed_time)
        data_from_dict_time_in_x.columns = [''] * len(data_from_dict_time_in_x.columns)
        return data_from_dict_time_in_x
        
    elif option_output == 'wave':
        #data_from_dict_wave_in_x = pd.DataFrame.from_dict(od, orient = 'columns')
        #data_from_dict_wave_in_x.index.name = 'relevant peaks'
        return data_from_dict_wave_in_x
        #print(data_from_dict_wave_in_x.to_string(), file = isf_results_fixed_wave)



#---------------------------------------------------------------------------------------------------------------------#
def positions_list_com_deleting_to_get_log_values(com_list2000, init_number, ending_number, step):
    
    returned_list = []
    def comparison_list(init_number, ending_number, step):
        comparison_list_x = []
        for i in range(init_number, ending_number, step):
            comparison_list_x.append(i)
        return comparison_list_x
    
    to_filter = comparison_list(init_number, ending_number, step)

    for row_in in com_list2000:
        #print(row_in)
        if int(row_in[0]) in to_filter:
            returned_list.append(row_in)
            
            
    return returned_list

def positions_list_com7(file_com):
    #THIS function exists because the intersection  of data in different files..
    coords = []
    x_l = []
    y_l = []
    z_l = []
    for skip in range(7): #Puse 7, solo xq el archivo empieza en 1000, en lugar de 2000
        next(file_com)
    for j, rows in enumerate(file_com):
        
        values = rows.split()
        if j%4 == 0:
            coords.append(int(values[0]))
        elif (j-5)%4 == 0:
            x_l.append(float(values[1]))
        elif (j-6)%4 == 0:
            y_l.append(float(values[1]))
        elif (j-7)%4 == 0:
            z_l.append(float(values[1]))
    
    xyz =   [list(x) for x in zip(x_l, y_l, z_l)]
    xyz_c = [list(m) for m in zip(coords, xyz)]

    return xyz_c

def get_log_values_for_no_com(isf_steps, step_to): #NOT USED HERE
    first_list, list_log_values_no_com = ([] for i in range(2))
    for i in isf_steps:
        first_list.append(int(i))
    
    for j in first_list:
        if j < step_to:
            list_log_values_no_com.append(j)
    #print(list_log_values_no_com)
    return list_log_values_no_com

def positions_list_com0(file_com):
    coords = []
    x_l = []
    y_l = []
    z_l = []
    for skip in range(0): #Puse 6, solo xq el archivo empieza en 1000, en lugar de 2000
        next(file_com)
    for j, rows in enumerate(file_com):
        
        values = rows.split()
        if j%4 == 0:
            coords.append(int(values[0]))
        elif (j-5)%4 == 0:
            x_l.append(float(values[1]))
        elif (j-6)%4 == 0:
            y_l.append(float(values[1]))
        elif (j-7)%4 == 0:
            z_l.append(float(values[1]))
    
    xyz =   [list(x) for x in zip(x_l, y_l, z_l)]
    xyz_c = [list(m) for m in zip(coords, xyz)]

    return xyz_c


def positions_list_com(file_com):
    coords = []
    x_l = []
    y_l = []
    z_l = []
    for skip in range(3): #Puse 6, solo xq el archivo empieza en 1000, en lugar de 2000
        next(file_com)
    for j, rows in enumerate(file_com):
        
        values = rows.split()
        if j%4 == 0:
            coords.append(int(values[0]))
        elif (j-5)%4 == 0:
            x_l.append(float(values[1]))
        elif (j-6)%4 == 0:
            y_l.append(float(values[1]))
        elif (j-7)%4 == 0:
            z_l.append(float(values[1]))
    
    xyz =   [list(x) for x in zip(x_l, y_l, z_l)]
    xyz_c = [list(m) for m in zip(coords, xyz)]

    return xyz_c
            
            
def positions_list(file_tx):
    dict_atoms = {}
    coords = []
    for skip in range(9):
        next(file_tx)
    for rows in file_tx:
        values = rows.split()
        keys = int(values[0])
        coords = [float(values[1]), float(values[2]), float(values[3])]
        dict_atoms[keys] = coords   
        
    positions = sorted(dict_atoms.items())    
    return positions


if __name__ == "__main__": main()
