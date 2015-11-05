# tau_fit
This is a new fit method for tau analysis, published in 2013:  
Evidence for the Appearance of Atmospheric Tau Neutrinos in Super-Kamiokande, Phys. Rev. Lett. 110, 181802  
  
The basic idea is to fit the tau normalization in the paper with all systematic errors simultaneously with RooFit.
The basic PDFs are built from simulation of e/mu events and tau events separately. Each systematic error has a
set of PDFs in the fit, which is based on the Osc3++.  
The framework of compiling RooFit code is copied from  
https://github.com/IPNL-CMS/MttTools  
The required libraries could be installed by compiling build-external.sh  

The fundamental PDFs of background and signal is built from NN output file. The PDFs for systematic errors is built
with fijTo2D.py.

plot_fit.py is added to make plot of fit result of systematic errors.
