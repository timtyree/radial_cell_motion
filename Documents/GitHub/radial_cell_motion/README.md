<h4>
* Python       :   Cellular Automata for Radial Chemotaxis Monte Carlo<br>
*<br>
* PROGRAMMER   :   Timothy Tyree<br>
* DATE         :   Fri 13 Dec 2019 <br>
* PLACE        :   Rappel Lab @ UCSD, CA<br>
</h4>



<p>This project simulates radial cell motion resulting from constant radial speed models of the form</p>

 v<sub>r</sub>( g<sub>A</sub>,  g<sub>R</sub>) = v<sub>0</sub> sign(a g<sub>A</sub> - b g<sub>R</sub>)

<p>where  g<sub>A</sub> and g<sub>R</sub> are gradients or fractional gradients in a radial concentration profile.
g<sub>R</sub> is taken for steady state for a constitutively produced chemorepellent, R, and g<sub>A</sub> is taken to be a chemoattractant such as a cyclic adenosine monophosphate (cAMP) radial concentration field varying in time according to values stored in a pandas.DataFrame object, df_camp.  sign is the sign function, which is taken to be zero for arguments suffiently close to zero.  g_thresh is the minimum gradient that can be sensed to produce directed cell motion.  The default value is previously suggested, and for a 15µm long Dictyostelium discoideum cell with a 0.1nM cAMP drop accross its length, g_thresh = 0.1/15.<br>
(Song, 2006: https://doi.org/10.1016/j.ejcb.2006.01.012)</p>

<p>This module contains the Cell() class, which represents cell motion predicted by the constant radial speed models. Cell() is a 1D cellular automaton for radial cell motion in the presence of a chemoattractant and a chemorepellent.</p>
	
<p>This module also supports running batches of simulations in parallel and measuring/storing values of interest for each batch. </p>

<p>I made this module for my physics PhD research on chemotaxis, and although it serves my purposes, it is reasonably well-documented and riddled with "#TODO: " comments describing simple and clear incremental improvements.  Test cases are reasonably well developed.  If you're interested in helping make this module more useful to more people, I strongly encourage you to fork me on GitHub at TimtheTyrant.  Have a nice day!</p>
