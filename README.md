<h1>
Simulations of Radial Cell Motion
</h1>

<p>This repository contains the code used to generate the figures and data associated with "Cell dispersal by localized degradation of a chemoattractant" (Karmakar, Tyree, Gomer, Rappel, 2021).  Inward/Outward cell motion was simulated using local cell rules that suppose constant cell speed according to</p>

![equation](https://latex.codecogs.com/gif.latex?v%28g_A%2C%20g_R%29%20%3D%20v_0%20%5Ctext%7Bsign%7D%28a%20g_A%20-%20b%20g_R%20%29)

<p>where  g<sub>A</sub> and g<sub>R</sub> are gradients or fractional gradients in a radial concentration profile.
g<sub>R</sub> is taken for steady state for a constitutively produced chemorepellent, R, and g<sub>A</sub> is taken to model the gradient of chemoattractant, cyclic adenosine monophosphate (cAMP).  We compared simulation results to experimental microscopy data of the social migration of Dictyostelium discoideum during starvation.  We found no significant evidence for the existance of the aforementioned chemorepellant, R.  However, our experimental results support cell dispersal mediated by the secretion of a protein species, phosophodiesterase, that locally degrades the chemoattractant, cAMP.</p> 
	
<p>The radial concentration fields vary in time according forward Euler integration of a reaction-diffusion equation.  sign is the sign function, which is taken to be zero for arguments suffiently close to zero.  g_thresh is the minimum gradient that can be sensed to produce directed cell motion.  The default value is previously suggested, and for a 15µm long Dictyostelium discoideum amoeba with a 0.1nM cAMP drop accross its length, g_thresh = 0.1nM/15µm ≈ 0.0067 nM/µm (Source: https://doi.org/10.1016/j.ejcb.2006.01.012)</p>

<p>This repository contains the Cell() class, which models cell motion predicted by the aforementioned constant radial speed models. Cell() is a 1D cellular automaton for radial cell motion in the presence of a chemoattractant and a chemorepellent.</p>
	
<p>This repository also supports running batches of simulations in parallel and measuring/storing values of interest for each batch. A framework is provided for gridlike searches of parameter spaces that utilizes an HTCondor backend.  The repository code is predominantly written in Python.</p>

<p>I feel the code in this repository is reasonably well-documented, however it is riddled with "#TODO: " comments describing simple and clear incremental improvements. Feel free to contribute!</p>
