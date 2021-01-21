<h4>
Simulations of Radial Cell Motion
</h4>

<p>This repository contains the code used to generate the figures and data associated with "Cell dispersal by localized degradation of a chemoattractant" (Karmakar, Tyree, Gomer, Rappel, 2021).  Inward/Outward cell motion was simulated using local cell rules that suppose constant cell speed according to</p>


- <img src="https://latex.codecogs.com/gif.latex?v(g_A, g_R) = v_0 \text{sign}(a g_A - b g_R )" />

<p>where  g<sub>A</sub> and g<sub>R</sub> are gradients or fractional gradients in a radial concentration profile.
g<sub>R</sub> is taken for steady state for a constitutively produced chemorepellent, R, and g<sub>A</sub> is taken to model the gradient of chemoattractant, cyclic adenosine monophosphate (cAMP) 
	
	The radial concentration field varyies in time according to values stored in a pandas.DataFrame object, df_camp.  sign is the sign function, which is taken to be zero for arguments suffiently close to zero.  g_thresh is the minimum gradient that can be sensed to produce directed cell motion.  The default value is previously suggested, and for a 15µm long Dictyostelium discoideum cell with a 0.1nM cAMP drop accross its length, g_thresh = 0.1/15.<br>
(Song, 2006: https://doi.org/10.1016/j.ejcb.2006.01.012)</p>

<p>This repository contains the Cell() class, which represents cell motion predicted by the constant radial speed models. Cell() is a 1D cellular automaton for radial cell motion in the presence of a chemoattractant and a chemorepellent.</p>
	
<p>This repository also supports running batches of simulations in parallel and measuring/storing values of interest for each batch. </p>

<p>I feel the code in this repository is reasonably well-documented, however it is riddled with "#TODO: " comments describing simple and clear incremental improvements. Feel free to contribute!  Have a nice day!</p>
