#!/usr/bin/env python3
# an overseer of the humble workers for the PDE Simulation.
# Developer: Tim (the Tyrant) Tyree
# 4.1.2020
# prepares data search grids and maneuvers data.  
# compresses output to a single .csv file

import pandas as pd
import numpy as np
from scipy.special import kn,kv #modified bessel function of the second kind of integer order n
from scipy.interpolate import BSpline

#automate the boring stuff
import time, os, sys, shutil
beep = lambda x: os.system("echo -n '\a';sleep 0.2;" * x)
if not 'nb_dir' in globals():
    nb_dir = os.getcwd()

##################################################
####Data handling functionality ##################
##################################################
def load_df(fn):
    '''load concetration simulation output.'''
    return pd.read_csv(fn)#, index_col='r')
def concatenate_results(df_1, df_2): 
    '''function that combines two custom pde_results .csv files, fn_1, fn_2, 
    and returns the resulting csv filename, fn_out.'''
    return pd.concat([df_1,df_2])

    #(deprecated) commented out is for adding time series data
    # times_1 = df_1.T.index.values.astype('uint16')
    # times_2 = df_2.T.index.values.astype('uint16')+times_1.max()
    # df_link = df_2.rename(columns=dict(zip(df_2.columns,times_2)))
    # df_out  = pd.concat([df_1.T,df_link.T]).T
    # return df_out
def get_last(df):
    '''extract the final field state from a given csv in fn.'''
    last_id = df.T.index.astype('uint16')[-1]
    df_out = df.T.loc[str(last_id)]
    return df_out
def print_args():
    print(sys.argv)
    return True

###########################################################
##go from input grid of inputs to dataframe of inputs######
###########################################################
def grid_to_df(
    kgrid = [0.02],
    Lgrid = [100.],
    c0grid= [4],
    Tgrid = [21],
    itnogrid = [2],
    dtgrid= [0.0005],
    trgrid= [10],
    modegrid = []):
    '''
    returns complete inputs as entries in a dataframe in correspondance to the input grids.
    specifying a list of input values for any of itnogrid, kgrid, 
    Tgrid, Lgrid, c0grid, trgrid will be iterated over.
    column_names = [ 'kPDE', 'LPDE', 'c0', 'T', 'iter_no', 'dt', 'time_res']
    ''' 
    column_names = [ 'kPDE', 'LPDE', 'c0', 'T', 'iter_no', 'dt', 'time_res', 'mode']
    sample_space = [[K ,L, c0, T, itno, dt, tr, mode]
                for mode in modegrid
                for K in kgrid
                for L in Lgrid 
                for c0 in c0grid
                for T in Tgrid
                for dt in dtgrid
                for itno in itnogrid
                for tr in trgrid
               ]

    sample_space = np.array(sample_space).transpose().astype('float16')
    return pd.DataFrame(dict(zip(column_names,sample_space)))

def df_to_dat(df,save_fn):
    '''df is the output of grid_to_df.
    save_fn is a string ending in .csv.'''
    df.to_csv(save_fn, header=False,sep=' ', index=False)
    shutil.copy(save_fn, save_fn.split('.')[0]+'.dat')
    os.remove(save_fn)



# desiered output format
# kPDE, LPDE, c0, T, iter_no=4, dt=0.0005, time_res = 10)


#TODO(later): make command line interface for for grid mesh generation
# if __name__ == "__main__":
# #  assign input parameter values
# #  assign defualt values from uniform random numbers
#     x_low = uniform(-10,0)
#     x_high = uniform(0,10)
#     y_low = uniform(-10,0)
#     y_high = uniform(0,10)
#     bound_array = [x_low, x_high, y_low, y_high]
#     for i in range(1,len(sys.argv)):
#         if i < 5:
# #  Replace the random values with the supplied values 
#             bound_array[i-1] = float(sys.argv[i])
   
# # The range for brute function requires in tuples  
#     brute_range = ((bound_array[0],bound_array[1]), (bound_array[2], bound_array[3]))
#     print('Search Boundary  x1= {0:3.3f} x2= {1:3.3f} x3= {2:3.3f} x4= {3:3.3f}'.format(*bound_array))
    
# # Here we are doing a brute force optimization. The function is evaluated in grids of points. 
# # brute_range is a tuple and defines the boundary for the grid points
# # finish=None means no local search. To make the search efficient choose finish=optimize.fmin
#     result_from_brute = optimize.brute(rosenbrock, brute_range, full_output=True, finish=None)
#     function_min = result_from_brute[1]
#     coordinate_of_min = result_from_brute[0]
#     #print ('Initial Coordinates= ',brute_range)
#     print ('Search Result= ',function_min, coordinate_of_min)

