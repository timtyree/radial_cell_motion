#!/usr/bin/env python3
# an overseer of the humble workers for the PDE Simulation.
# Developer: Tim (the Tyrant) Tyree
# 4.1.2020
# prepares data search grids and maneuvers data.  
# compresses output to a single .csv file
import os, sys, argparse, pandas as pd
# , numpy as np

##########################################
### Functionality for output collection ##
##########################################
def out_to_df(output_fn):
    '''takes a single job.out file and returns a df with the results.'''
    out_fn = output_fn
    os.path.exists(out_fn)
    h = open(out_fn)
    line = h.readline()
    #the first line should read "Printing Inputs".  If it doesn't, reject the whole file.
    if not ('Printing Inputs:\n'==line):
        return None
    #read the input line
    line = h.readline()
    str_list = line[:-1].split(' ')
    input_dict = {s.split(':')[0]:float(s.split(':')[1]) for s in str_list}
    #next couple lines not needed
    line = h.readline()
    line = h.readline()
    assert(line=='Printing Outputs including mean cell direction (mcd):\n')
    #read the output lines as dataframe
    lines = h.readlines()
    tmp_fn = out_fn.split('.')[0]+'_tmp.csv'
    open(tmp_fn, 'x').close()
    w = open(tmp_fn, 'w')
    w.writelines(lines)
    w.close()
    df = pd.read_csv(tmp_fn)
    os.remove(tmp_fn)
    #add inputs to dataframe
    for key,value in zip(input_dict.keys(),input_dict.values()):
        df[key] = value
    return df

def condense_output(save_fn = 'results.csv', output_dir = 'Log/'):
    os.chdir(output_dir)
    #get list of output files
    out_fn_list = []
    for fn in os.listdir():
        if (fn.split('.')[1]=='out'):
            out_fn_list.append(fn)      
    #make a single dataframe for all values in out_fn_list
    df_out = None
    for out_fn in out_fn_list:
        df = out_to_df(out_fn)
        df_out = pd.concat([df_out,df])
    df_out.reset_index(inplace=True)
    os.chdir('..')
    df_out.to_csv(save_fn)
    return True

#############################
### command line interface (that's as simple possible)
#############################
if __name__ == "__main__":
    # Construct the argument parser
    ap = argparse.ArgumentParser()
    # Add the arguments to the parser
    ap.add_argument("-o", "--foperand", required=True,
       help="output directory")
    ap.add_argument("-s", "--soperand", required=True,
       help="save filename ending in .csv")
    args = vars(ap.parse_args())
    output_dir = args['foperand']
    save_fn = args['soperand']

    #collect jobs
    print(f'saving the job.out files to the single file, {save_fn}.')
    condense_output(save_fn = save_fn, output_dir = output_dir)
    print("saving complete.")