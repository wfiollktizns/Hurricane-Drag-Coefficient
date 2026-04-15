# Hurricane-Drag-Coefficient
Matlab code for calculation of drag coefficient and related variables using saildrone data

Main program is saildrone_direct_flux_clean.m. Other m-files are required to run it.
Data required are netcdf files of the form sdnumTCNAMEres.nc, where 'sdnum' is a 4 or 5 digit saildrone number, 'TCNAME' is a tropical cyclone name (e.g., 'SAM'), and 'res' is the temporal resolution (1min,20hz,4hz). Drag coefficient and related variables are written to a mat-file within saildrone_direct_flux_clean.m.
