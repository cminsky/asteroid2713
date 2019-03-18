To reproduce the analysis, run these notebooks in order:
1. [DataReduction.ipynb](DataReduction.ipynb): This is where I calibrated my images, correcting for readout signal, dark current, and imperfections in the optical path.
2. [DataProcessing.ipynb](DataProcessing.ipynb): This is where I turned my calibrated images into lists of the asteroid's apparent magnitudes, using the classes defined in [Analysis.py](Analyis.py) which are for storing and analyzing FITS files and sets of FITS files from a single night of observation pointing at the same field of view.
3. [LightCurve.ipynb](LightCurve.ipynb): This is where I created the asteroid light curve, tried different rotational periods, and calculated the asteroid's primary axis ratio from the light curve amplitude.

All data were taken at MIT Wallace Astrophysical Observatory. I took the bias, dark, and light frames. The median flat was created by Dr. Amanda Bosh. The MIT 12.410 F18 teaching staff (Richard Binzel, Amanda Bosh, Michael Person, Jane Connor, and Ani Chiti) taught me how to phase a light curve and use photutils for photometry.
