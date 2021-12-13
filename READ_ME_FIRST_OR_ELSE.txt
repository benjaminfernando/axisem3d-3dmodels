This folder contains a series of important things. These are your standing orders. Read them and understand them!! 

1. The source code: this is just a copy of that at https://github.com/kuangdai/AxiSEM-3D/tree/master/SOLVER/src 
2. A template input file. This contains: 

- An example HPC submission script 
- A set of input files configured for a run with both surface and Moho topography, with GSN and YKR stations, and an earthquake in Guatemala. The exemplar 3D model files: 
	- moho.nc: 1 degree sampling, smoothed using a Gaussian with 1 degree standard deviation. The selected input file is not smoothed at the poles but you can easily add this in if you want to 
	- etopo01.nc: 1 degree sampling, also smoothed using a Gaussian with 1 degree standard deviation 

3. The codes* needed to create the 3D model input files for simulations like these, as used to create the 3D model files in the template_input directory
4. The raw data files needed to create these 3D model files, i.e. the raw copies of ETOPO1 and Crust 1.0's Moho
5. A sort of manual
6. A python script* to read in synthetic data from AxiSEM3D runs and plot them.

* Note that these python scripts are not designed to be elegant/pythonic, but rather easy to understand and edit. They are perhaps overly-verbosely commented, and in the 3D file scripts (moho.py and etop01.py) they include functionality that is perhaps not needed (e.g. smoothing at the poles), but was something that I wrote to try and solve stability issues. You can just comment it out, set the relevant flags to false, or delete it; but I thought I would include it just in case you want to build upon any of this. 
