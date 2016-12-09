# spiking_variability
Simulation codes for https://arxiv.org/abs/1605.08909

To reproduce figures, have the tools folder on your MATLAB path. Tested on OS X and Linux. 

**Figures 8 and 10**
```
gen_figs_8_10(50,1.5,6,'seed',2);
```
to replot figures once the simulation has been run:
```
gen_figs_8_10(50,1.5,6,'seed',2,'loadSavedSim',true);
```

**Figure 9**
```
gen_fig_9(true)
```
to replot figures once the simulation has been run:
```
gen_fig_9(false)
```
