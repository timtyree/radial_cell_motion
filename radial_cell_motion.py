# /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  * Python       :   Cellular Automata for Radial Chemotaxis Monte Carlo
#  *
#  * PROGRAMMER   :   Timothy Tyree
#  * DATE         :   Fri 13 Dec 2019 
#  * PLACE        :   Rappel Lab @ UCSD, CA
#  *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  */

'''
This project simulates radial cell motion resulting from constant radial speed models of the form

vr = v0*sign(a*ga-b*gr), 

where ga and gr are gradients or fractional gradients in a radial concentration profile.
gr is taken to be at steady state for a constitutively produced chemorepellent, R, and ga is taken to be a chemorepellent
such as a cyclic adenosine monophosphate (cAMP) radial concentration field varying in time according to values stored in
a pandas.DataFrame object, df_camp.
'''

import pandas as pd, numpy as np
from scipy.special import kv
from scipy.interpolate import BSpline
from joblib import Parallel, delayed
import multiprocessing
num_cores = multiprocessing.cpu_count()

#automate the boring stuff
import time, warnings, sys, os

#TODO: make automated database repository creation function for results to stream to
#TODO: stream simulated results to that database repository as they're made so they don't take up virtual memory
#TODO: implement method handling for 'linlin', 'linlog', 'loglin', and 'loglog'
#TODO: make (an np.vectorized) get_cell_trajectory function
#TODO: implement method kwrd in simulate
#TODO: incorporate minimum gradient sensing into Cell
#TODO: repeat for pde degradation model outward gradients
#TODO: put test cases in their own file and have an associated test df_camp to pull without user input df_camp

# import the library
transfer_dir = os.path.join('pde_sim/nb')
if transfer_dir not in sys.path:
    sys.path.append(transfer_dir)
transfer_dir = os.path.join('pde_sim/pde-sim-transfer')
# if transfer_dir not in sys.path:
#     sys.path.append(transfer_dir)
# import pde_sim_returns_camp as pdesim
# # import pde_sim.nb.pde_sim_returns_camp as pdesim


############################################################################
# Cell Class and helper functions
############################################################################

def interpolate(x,y,n=3):
	'''return get order n spline interpolation function fit to array x and array y.'''
	return BSpline(x, y, n, extrapolate=True)

def sign_w_thresh(arg,thresh):
	    '''returns 1 if arg is more than thresh, -1 if arg is less than -thresh, and zero else.
	    assumes thresh>0. See test cases for example usage.'''
	    if arg-thresh>0:
	        return 1
	    elif arg+thresh<0:
	        return -1
	    else:
	        return 0

class Cell():
	'''1D cellular automaton for radial cell motion in the presence of a chemoattractant and a chemorepellent.
	vr = v0*np.sign_w_thresh(a*GA-b*GR), where GA, GR are some gradient functions.
	radial unit vector is pointed away from the origin.
	All these cell properties might slow things down.  Have caution.
	g_thresh is the minimum gradient that can be sensed to produce directed cell motion.  The default value 
	is previously suggested as a 15µm long Dicty. cell in a 0.1nM cAMP drop accross its length. 
	(^Song, 2006: https://doi.org/10.1016/j.ejcb.2006.01.012)
	TODO: incorporate cell memory?
	'''
	def __init__(self, r=50, a=1, b=1, v0=5, dt=1./6., rmin=50,rmax=1000,time=0, record_traj=True, g_thresh=0.1/15):
		self.v0 = v0 #velocity scale (µm/min)
		self.a  = a #sensitivity to cAMP
		self.b  = b #sensitivity to R
		self.r  = r #current location µm
		self.start_value = r #initial location µm
		self.dt = dt#time step in minutes
		self.time     = time #last recorded time
		self.in_rmeshQ= self.in_rmesh(rmin,rmax)
		self.g_thresh = g_thresh#nM/µm min gradient that can be sensed to produce directed motion.
		self.traj_lst = [r]
		return None
	
	def get_r(self):
		return self.r
	
	def set_r(self, r):
		self.r = r
		return self
	
	def calc_vr(self, GA, GR):
		'''return the velocity given by v0*np.sign(a*GA - b&GR).
		GA and GR is the local value of chemoattractant and chemorepellant concentration gradient in nM/micron'''
		return self.v0*sign_w_thresh(self.a*GA - self.b*GR,self.g_thresh)
		# return self.v0*np.sign(self.a*GA - self.b*GR)
		
	def move(self,GA,GR):
		'''radially move by an amount vr*dt'''
		self.r = self.r + self.calc_vr(GA,GR)*self.dt
		self.traj_lst.append(self.r)
		self.time = self.time + self.dt
		return self
	
	def in_rmesh(self,rmin,rmax):
		'''checks if cell is still in rmesh.'''
		boo = (rmin<=self.get_r()) & (rmax>=self.get_r())
		if not boo:
			self.in_rmeshQ = boo
		return boo


############################################################################
# Routine for Simulating Cell Motion
############################################################################

def simulate(Rmax,LR,b, df_camp, method='linlin', g_thresh=0.1/15):
	'''Returns list of Cell objects after integrating in time.  
	Simulates cell motion according to a binary model for radial cell motion.
	Rmax = max concentration of chemorepellent in nM
	LR   = characteristic lengthscale of chemorepellent in µm
	b    = relative sensativity of chemorepellent to chemoattractant.
	df_camp   = a pandas.DataFrame with radial camp profiles progressing through time
		'r' column depicting an equally spaced mesh of distances from the origin that has the form
		many more columns with int headers representing the time (in seconds) since start of camp simulation
		for example, df_camp[10] = the radial camp concentration 10 seconds after the start of simulation.
	Set Rmax to zero to ignore the effects of any chemorellent.
	'''
	df = df_camp
	r_values = df['r'].values
	dt   = np.diff(df.columns[1:].values).mean()/60#time step in minutes
	r0   = r_values[0]
	rmin = r0
	rmax = r_values[-1]
	method   = method.lower()



	# calculate steady state chemorepellent field and interpolate
	if Rmax >= 0:
		#use either gradient or fractional gradient of constitutivley produced chemorepellent
		if (method=='linlin') | (method=='loglin'):
			#calculate R gradient field
			DR = -Rmax/LR*kv(1,r_values/LR)/kv(0,r0/LR)
			dr = interpolate(r_values,DR)
			gr = dr
		elif (method=='linlog') | (method=='loglog'):
			#calculate R fractional gradient field
			with warnings.catch_warnings(record=True) as w:
				DlnR= -1/LR*kv(1,r_values/LR)/kv(0,r_values/LR)
				#replace dividebyzero warnings with an asymptotic fractional gradient 
				DlnR[np.isnan(DlnR)] = -(1+LR/2/r_values[np.isnan(DlnR)])/LR 
				#turn of DlnR if there Rmax is zero.
			dlnr = interpolate(r_values,DlnR)
			gr   = dlnr
		else:
			pass
	else: 
		DR = 0*r_values
		dr = interpolate(r_values,DR)
		gr = dr
	
	if (method not in ['linlin','linlog','loglin','loglog']):
		warnings.warn('''"Unexpected Method at RLBM = ({},{},{},{}): 
			using zero for chemorepellent and linear gradient 
			for chemorepellent."'''.format(Rmax,LR,b,method))
		DR = 0*r_values
		dr = interpolate(r_values,DR)
		gr = dr	
		method = 'linlin'
	
	#integrate and view the trajectory of a multiple cells
	start_values = np.arange(r0+5,800,5)
	cell_lst = []
	for r in start_values:
		cell_lst.append(Cell(r=r, dt=dt, b=b, g_thresh=g_thresh))

	#for each timestep
	inputs   = df.columns[1:].values
	for tQ in inputs:

		# calculate current chemoattractant field and interpolate
		if (method=='linlin') | (method=='linlog'):
			#calculate A gradient field
			DA = df[tQ].diff().values
			da = interpolate(r_values,DA)
			ga = da
		elif (method=='loglin') | (method=='loglog'):
			#calculate A fractional gradient field
			with warnings.catch_warnings(record=True) as w:            
				DlnA = df[tQ].diff().values/df[tQ].values
				#replace dividebyzero warnings with zero
				DlnA[np.isnan(DlnA)] = 0*r_values[np.isnan(DlnA)]
				dlna = interpolate(r_values,DlnA)
				ga   = dlna

		#for each cell still in r_values
		for cell in cell_lst:
			if cell.in_rmeshQ:
				#move and record motion
				r = cell.get_r()
				cell.move(GA=ga(r), GR=gr(r))
				#if a cell leaves rmesh, update cell.in_rmeshQ
				cell.in_rmesh(rmin,rmax)
			else:
				pass
	return cell_lst


############################################################################
# Routine for Measuring/Summarizing Cell Motion
############################################################################

#measure if any cells are attracted, and if so, what is the max range of net attraction
@np.vectorize
def increasedQ(cell):
	'''returns true if the final value of cell.traj_lst is greater than the initial value of the array.
	arguement should be similar to cell_lst[0].traj_lst.'''
	return cell.r>cell.start_value

@np.vectorize
def turnedQ(cell):
	'''returns true if the cell changed directions at all during its journey.'''
	traj = np.array(cell.traj_lst)[~np.isnan(cell.traj_lst)]
	return np.diff(np.diff(traj)>0).any()

@np.vectorize
def get_persistance_time(cell):
	'''measure inward persistance time for a given cell.'''
	if np.isnan(cell.traj_lst).any():
		#instead of removing any np.nan values, just return np.nan
		#traj = np.array(cell.traj_lst)[~np.isnan(cell.traj_lst)]
		return np.nan
	else:
		moving_outward = np.diff(cell.traj_lst)>0
		turned_inward  = moving_outward.argmin()
		turned_outward = len(moving_outward) - moving_outward[::-1].argmin()
		time_spent     = (turned_outward-turned_inward)*cell.dt
		return time_spent

def measure(cell_lst):
	'''returns a dict of values summarizing the results in cell_lst.
	max_range is the furthest distance where cells dispersed overall, and is np.nan if no dispersal occurred.
	max_range_felt is the furthest distance where cells changed direction at all, and is np.nan if no change of direction occurred.
	pt_mean is the mean persistance time of cells that turned inward and then turned outward.
	I assume here that at most one cAMP pulse occurs and that each cell turns at most twice.
	'''
	boo = increasedQ(cell_lst)
	any_dispersed = bool(boo.any())
	if any_dispersed:
		id_furthest   = len(boo)-boo[::-1].argmax()-1
		max_range     = int(cell_lst[id_furthest].get_r())
	else:
		max_range     = np.nan 

	boo = turnedQ(cell_lst)
	any_turned = bool(boo.any())
	if any_turned:
		id_furthest   = len(boo)-boo[::-1].argmax()-1
		max_range_felt= int(cell_lst[id_furthest].get_r())
		#measure inward persistance time for all cells
		persistance_times = get_persistance_time(cell_lst)
		#summarize inward persistence times
		pt = np.array(persistance_times)[~np.isnan(persistance_times) & boo]
		pt_min, pt_max, pt_mean = (float(np.around(pt.min(),2)),float(np.around(pt.max(),2)),float(np.around(pt.mean(),2)))
	else:
		max_range_felt= np.nan 
		pt_min, pt_max, pt_mean = (np.nan,np.nan,np.nan)
	output_trial = [any_dispersed, max_range, any_turned, max_range_felt, pt_min, pt_max, pt_mean]
	columns_trial = ['any_dispersed', 'max_range', 'any_turned', 'max_range_felt', 'pt_min', 'pt_max', 'pt_mean']
	out_dict = dict(zip(columns_trial,output_trial))
	return out_dict



############################################################################
# Batch Simulator
############################################################################

def run_batch(Rmaxgrid, LRgrid, bgrid, df_camp, method_lst=['linlin'], g_thresh=0.1/15):
	'''run a batch over the input grids in parallel.
	method is a list of methods to be iterated over.
	valid methods are "linlin", "linlog", "loglin", and "loglog".'''
	#check for user error
	for method in set(method_lst):
		if (method.lower() not in ['linlin','linlog','loglin','loglog']):
			print("Unexpected Method: "+method, sys.exc_info()[0])
			raise
	if type(df_camp)!=pd.DataFrame:
		print("df_camp has type {} instead of pandas.DataFrame.".format(type(df_camp)))
		raise
	
	#build sample_space and process inputs
	sample_space = [(Rmax,L,b,method) for Rmax in Rmaxgrid for L in LRgrid for b in bgrid for method in tuple(set(method_lst))]
	inputs = range(len(sample_space))
	def processInput(item):
		columns_batch = ['Rmax', 'LR','b', 'method']
		Rmax   = sample_space[item][0]
		LR     = sample_space[item][1]
		b      = sample_space[item][2]
		method = sample_space[item][3].lower()
		try:
			cell_lst    = simulate(Rmax,LR,b, df_camp=df_camp, method=method, g_thresh=g_thresh)
			batch_dict  = dict(zip(columns_batch,[Rmax,LR,b,method]))
			output_dict = measure(cell_lst)
			for K in batch_dict.keys():
				output_dict[K] = batch_dict[K]
			return output_dict
		except Exception as e: 
			return e, 'Error occured at input item RLBM = ({},{},{},{}).'.format(Rmax,LR,b,method)

	start = time.time()
	print('simulation starting...')
	lst_output = Parallel(n_jobs=num_cores, verbose=.2)(delayed(processInput)(item) for item in inputs)
	end   = time.time()
	print('{} seconds elapsed running simulation.'.format(np.around(end-start,1)))
	return lst_output

############################################################################
# Test Cases
############################################################################

def run_test_cases():
	'''test cases for the contents of radial_cell_motion.py.
	Some of these are redundant, but in a way that aids in error identification.'''
	#load test DataFrame
	test_data_dir = 'data/mg_source_spatial_propagation_pde_deg_itno_2_K_0_T_6_L_60_c0_0.0_10s_obs.csv'
	df_camp = pd.read_csv(test_data_dir)
	#test cases for sign function with thresholding
	assert(sign_w_thresh(1--1,1.9)==1)
	assert(sign_w_thresh(1--1,2)==0)
	assert(sign_w_thresh(-1-1,0)==-1)
	assert(sign_w_thresh(np.nan,999)==0)
	#test cases for Cell
	assert(Cell())
	cell = Cell(r=100,dt=1)
	assert(cell.move(-1,0).get_r()==95)#move towards cAMP
	assert(cell.move(0,-1).move(0,-1).get_r()==105)#move away from R
	assert(cell.move(0,0).get_r()==105)#no gradients mean stay still
	assert(Cell().dt*6==1)#assert time step inititialized right (10seconds)
	assert(cell.in_rmesh(100,110))
	assert(~cell.in_rmesh(100,90))
	assert(~cell.in_rmesh(150,110))
	#test cases for simulation
	assert(simulate(Rmax=-1,LR=1,b=1,df_camp=df_camp))
	assert(simulate(Rmax=0,LR=1,b=1,df_camp=df_camp))
	assert(simulate(Rmax=1,LR=1,b=1,df_camp=df_camp))
	# assert(simulate(Rmax=1,LR=1e-8,b=1))#just don't do small LR for now
	cell_lst = simulate(Rmax=1,LR=1,b=1,df_camp=df_camp);
	out_dict = measure(cell_lst)
	#TODO: test cases for measure
	cell_lst = simulate(Rmax=1,LR=1,b=1,df_camp=df_camp);
	out_dict = measure(cell_lst)
	#test cases for batch simulator
	batch_dict = run_batch([1],[1],[1,2,3],df_camp=df_camp);
	assert(batch_dict!=None)
	assert(type(batch_dict)==list)
	assert(len(batch_dict)==3)
	assert(type(batch_dict[0])==dict)
	assert(len(batch_dict[0].keys())>=10)
	assert(type(batch_dict[0]['any_dispersed'])==bool)
	assert(type(batch_dict[0]['any_turned'])==bool)
	#TODO: add type test cases for nonboolean output of measure
	#test different methods
	assert(type(measure(simulate(1,1,1,df_camp,'linlin')))==dict)
	assert(type(measure(simulate(1,1,1,df_camp,'linlog')))==dict)
	assert(type(measure(simulate(1,1,1,df_camp,'loglin')))==dict)
	assert(type(measure(simulate(1,1,1,df_camp,'loglog')))==dict)
	#test method batching
	df_out = pd.DataFrame(run_batch([1],[1],[1],df_camp,method_lst=['loglog', 'Linlin', 'LOGLIN', 'linLoG', 'loglog']));
	assert(len(list(df_out['method'].values))==4)
	method_lst = {'linlin', 'loglog', 'linlog', 'loglin'}
	assert(set(list(df_out['method'].values))==method_lst)
	#verify different g_thresh batches give different results
	tmp1 = run_batch([1], [1], [1], df_camp = df_camp, method_lst=method_lst, g_thresh=0);
	tmp2 = run_batch([1], [1], [1], df_camp = df_camp, method_lst=method_lst, g_thresh=1000);
	assert(~pd.DataFrame(tmp1).equals(pd.DataFrame(tmp2)))
	#TODO: make run_batch robust to repeated entries in method_lst
	print('all test cases passed')




