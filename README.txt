This is the readme for the models associated with the paper:

Handy G, Taheri M, White JA, Borisyuk A (2017) Mathematical
investigation of IP3-dependent calcium dynamics in astrocytes.
J Comput Neurosci

These MATLAB and XPP files were contributed by G Handy.

Use example_ca_sim_TH.m to reproduce calcium simulations found in
Fig. 2 

Use autoRespTypeDetection_TH.m to generate calcium responses 
resulting from multiple (or all 600) underlying IP3 traces, as 
well as characterize those calcium and IP3 traces (in terms of, 
for example, amplitude, total duration, area under the curves, etc.).

Use ca_bifurcation_xpp.ode with XPPAUT to reproduce the bifurcation
diagrams found J. Computation Neuroscience (submitted)

Supporting_Functions contains all necessary functions to run the
simulations.

Histogram_Figures contains the data base of calcium response types
distributions for a range of parameter values (original figures can be
found in J. Computation Neuroscience (submitted))

Hopefully you have found this code helpful!

--Marsa Taheri and Gregory Handy
