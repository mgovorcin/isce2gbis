# isce2gbis
Function to generate GBIS INSAR input file from ISCE generated interferogram

 Usage: GBISinput = ISCE2GBIS(isce_dir, seed_ref, crop)
 isce_dir  = '/dir/' - define the directory with filt_topophase.unw.geo file

 output_dir = 'dir' -define the output dir for GBIS insar input file

 seed_ref = [Lon,Lat] - define the unwrapping seed point (default left upper corner), move it closer to deformation area but not on the area affected by deformation (OPTIONAL)

 crop = [lon_min, lon_max, lat_min, lat_max] - Define the window for
 
 # cropping (OPTIONAL)
 =========================================================================
 This function is created by: Marin Govorcin (22.08.2019)
 
 Address: Faculty of Geodesy, University of Zagreb
 
 Email: mgovorcin@geof.hr
 
 TO DO LIST: Script works with StripmapApp.py- update it to work with TopsApp.py and InsarApp.py (code lines: 99, 102, 138)
