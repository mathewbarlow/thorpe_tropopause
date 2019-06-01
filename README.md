# thorpe_tropopause

<b> Overview </b>

This repository contains preliminary python code to replicate the results of Thorpe (1985), Thorpe (1986), and Fig. 15a from Hoskins et al. (1985), which show the circulation for various axisymmetric distributions of potential vorticity and boundary variations in potential temperature.  There is a separate program for each figure considered, although the only changes are in the given conditions and parameter values.

Appears to work for all the figures from Thorpe (1985) EXCEPT Fig. 3, the sloping tropopause case, although the results are very similar if the program is halted before convergence (as determined in terms of potential vorticity). Does not work yet for the figures from the other papers, although very similar results can be obtained by doubling the radial scale of the given conditions and halting before convergence.

Comparisons between the output of the programs and the original figures are provided in the "figures" folder of this repository.

<b> Summary of Results </b>

The code appears to work fine for Figs. 1, 2, 4, 5, and 6 from Thorpe (1985).  For example, here is a comparison for Fig. 5, which shows a close match. (Note that there appears to be a typo in the paper, it looks like the sine should be squared in the specification of the potential vorticity, which is what done here.)

<img src="figures/comparison_thorpe_1985_fig_5.png" width="500" height="500">

<b> Discussion </b>

For simplicity, potential vorticity and potential temperature are linearly interpolated when determining the tropopause.  I have tried more sophisticated approaches but they didn't seem to make much difference and tended to cause convergence problems in some cases, so I did not include them here.  I also did some sensitivity tests with much higher resolution, and without any interpolation for potential vorticity (every grid box either has the tropospheric or stratospheric value with no intermediate values), and the way the tropopause is treated in the code does not appear to explain the differences from the original figures.  

Given the similarity between the results and the original figures when I double the radial scale of the given conditions, I've either misunderstood something about the basic parameters, how to implement them, or have a bug somewhere. I've traded various permuations of different values but can't figure out the problem.

<b> Acknowledgements </b>

Support from NSF AGS-1623912 and NSF AGS-1657921 is gratefully acknowledged.

<b> References </b>

Hoskins, B.J., McIntyre, M.E. and Robertson, A.W., 1985. On the use and significance of isentropic potential vorticity maps. Quart. J. Roy. Meteor. Soc., 111(470), pp.877-946.

Thorpe, A., 1985: Diagnosis of balanced vortex structure using potential vorticity.  J. Atmos. Sci., 42, 397-406.
