# Neural-AP-morphofilt
Matlab files for simulating extracellular action potentials by morphological filtering (see paper, draft submitted to JCN). 
Python/LFPy files are also provided for comparision purposes.

These files reproduce figure 8 from the paper by launching demo_github.m (tested in Matlab 2018a). By changing the parameters 
in the script (axon and dendrite lengths and diameters, as well as the propagation speed, somatic dipole weigth and orientation), 
other bal-stick type morphologies can be reproduced as well. 
The main function is morphofiltd.m (hhrun.m is an implementation of the HH model).
The folder is self contained, you just need to put it on your Matlab path, or to browse to the JCN_demo/Matlab folder.

A Python folder is also provided, with a NEURON .hoc file inside. In order to launch the demo_github.py file, you need to have 
LFPy and NEURON installed (it was tested with versions 2.0.2 and 7.6.4 respectively), Python 3.5.
You can change the morphology also by loading a different .hoc.
