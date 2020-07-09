#!/usr/bin/env python3
import numpy as np, pandas as pd, os, sys, matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import


plt.rc('axes', linewidth=2)
#automate the boring stuff
import time, os, sys
beep = lambda x: os.system("echo -n '\a';sleep 0.2;" * x)
if not 'nb_dir' in globals():
    nb_dir = os.getcwd()
print('notebook is at: ' + nb_dir)
#TODO: make the abstraction
#DONE: make it greyscale
#TODO: try s = c0 and c = any output
def plot_for(df_raw,constant,x_col,y_col,c_col,
    s_col=None, s_scale=200,
    x_label = '$\mathbf{k_{PDE}}$ (1/s)',
    y_label = '$\mathbf{L_{PDE}}$ (µm)',
    c_label = '$\mathbf{C_{back}}$ (nM)',
    save_dir = f'{nb_dir}/../fig/', 
    y_log_scale=True, 
    x_log_scale=True, 
    axis_limits = [0.0001, 1, 10, 1000],
    fontsize = 18,
    show_cbar = True,
    marker = 's',
    cmap = None,
    alpha = 1,
    save_folder = f'{nb_dir}/../fig/',
    fig_size = [6.2, 3.5],
    c_limits = None,
    edgecolors = 'face',
    **kwargs):
    '''plot in 2D for constant-variable-response from grid search results for an abstract experiment.
    kwargs passed to matplotlib.pyplot.scatter(**kwargs).'''
    #step 1: select/compute the data to be plotted
    cond = f'mean_c<1000'#default conditions
    save_fn = f'plot_for constant'
    for x in list(constant):
        cond = cond + f' and {x}=={constant[x]}'#for each user defined condition
        save_fn = save_fn + f' {x}_{constant[x]}'
    df = df_raw.query(cond).copy()#[['kPDE','LPDE', 'c0', 'mean_c', 'dispersedQ']]
    save_fn = save_fn + f' variable_{x_col}_{y_col} response_{c_col}.png'
    
    # target_lst = ['mcd_at_100']#['mean_c']#['dispersedQ']#
    # c_lst = ['gray']
    # s_lst = [1]#zorder
    # fill_above_lst = [False]
    
    os.chdir(save_folder)
    if not cmap:
        cmap = plt.cm.get_cmap('gray_r')

    #step 2: plot the data
    #plot log scale phase diagram 
    fig = plt.figure(figsize=fig_size)#figure size
    ax  = fig.add_subplot(111)

    ##for each target/category 
    # for target, c, s, fill_above in zip(target_lst, c_lst, s_lst, fill_above_lst):

    #compute the features 
    C = df[c_col].values
    X = df[[x_col, y_col]].values
    if not s_col:
        S = s_scale
    else:
        S = s_scale*df[s_col].values
    im = ax.scatter(x=X[:,0],y=X[:,1], c=C, s=S,
    marker = marker, cmap=cmap, alpha = alpha, edgecolors = edgecolors)#plt.cm.jet)#
    # im = ax.scatter(x=X[:,0],y=X[:,1], c=C, s=S,
    # marker = marker, cmap=cmap, alpha = alpha)#plt.cm.jet)#
    # # set color limits as follows,
    if not c_limits:
        im.set_clim(df[c_col].values.min(), df[c_col].values.max())
    else:
        im.set_clim(*c_limits)

    if show_cbar:
        #color bar
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)#, labelsize=fontsize)
        cbar.set_label(c_label, fontsize=fontsize)#, fontweight='bold')
        cbar.ax.tick_params(width=2, length=6, which='both', direction='inout', labelsize=fontsize-4)

    #format plot
    # ax.set_title('background cAMP initially = {} nM'.format(c0), fontsize=fontsize, fontweight='bold')
    ax.set_xlabel(x_label, fontsize=fontsize, fontweight='bold')
    ax.set_ylabel(y_label, fontsize=fontsize, fontweight='bold')   
    
    if x_log_scale:
        ax.set_xscale('log')
    if y_log_scale:
        ax.set_yscale('log')
    if axis_limits is not None:
        ax.axis(axis_limits)
    # make ticks big and bold
    ax.tick_params(width=2, length=6, which='both', direction='inout', labelsize=fontsize)
    for tick in ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(fontsize)
        tick.label1.set_fontweight('bold')
    os.chdir(save_dir)    
    fig.savefig(save_fn, bbox_inches ='tight', dpi = 600)
    return fig