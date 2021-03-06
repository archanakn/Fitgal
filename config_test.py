"""Module to set fitgal config options and parameters"""
import numpy as np
from astropy.io import fits
import pylab
import numpy.ma as ma
# Select which components to fit
fit_bulge=True
fit_disk=True
fit_point=False
do_convolve=True
niter=0
#nxpts
#nypts
# Now set the initial value,mixnvalue,maxvalue,stepsize,isfixed for parameters
# isfixed=1 implies parameter is fixed, set isfixed=0 for non-fixed parameters
i0bvals = (2500.0,300,5000.0,10.1,0)
revals= (6.,1.0,30,0.1,0)
ebvals= (0.3,0.1,0.5,0.01,0)
nvals=  (6.0,1.0,6,0.01,0)
i0dvals= (1500,500,2000.0,10.,0)
rdvals= (10.0,1.0,10.0,0.1,0)
edvals= (0.8,0.1,0.8,0.01,0)
pointvals= (0.0,-1000.0,10000.0,10.0,1)
bpavals= (35.0,25.0,40.0,0.1,0)
dpavals= (35.0,25.0,40.0,0.1,0)
backgroundvals= (0.0,-100.0,2000.0,0.1,0)


#i0bvals = (4700.0,1000,50000.0,0.1,0)
#revals= (35.0,1.0,70,0.1,0)
#ebvals= (0.25,0.0,0.8,0.01,0)
#nvals=  (4.0,1.0,6,0.1,0)
#i0dvals= (200,0.0,1000.0,1.0,0)
#rdvals= (4.0,0.1,10.0,0.1,0)
#edvals= (0.3,0.0,0.8,0.01,0)
#pointvals= (0.0,-1000.0,10000.0,10.0,1)
#bpavals= (71.0,10.0,75.0,0.1,0)
#dpavals= (71.0,10.0,75.0,0.1,0)
#backgroundvals= (0.0,-100.0,2000.0,0.1,0)


have_pylab=True
# Do we have pylab? If yes, variable is set to True automatically
# by the program. It is best to leave this unchanged.

make_plots=True
# Set make_plots=False if you want to disable plot generation.

use_mask=True
# Set to True if you need to use mask files
# Currently, the mask file must be the same size as the image
# its value should be zero for valid pixels and 1 for invalid pixels

# WWW directory where fitgal.html is being generated
webdir='http://meghnad.iucaa.ernet.in/~yogesh/fitgal/' # trailing slash required!

# Parameters that control what outputs are generated and what they look like
plotiter = 10 # Number of iterations after which progress plot will be written
make_movie=True # Set True if you want to generate an animated GIF of fitting

# The following parameters should normally be provided as command line
# arguments. If the right number of command line arguments are not
# provided by the user the following default values will be used.

#Set the xcenter and ycenter of the galaxy in pixel units
xcenter=26.74 # test1
ycenter=27.53







#xcenter= 27.28 # test2
#ycenter= 26.59

#xcenter= 25.43 # test3
#ycenter= 24.29

# List the input data files to use
imagefile = 'model_galaxy.fits' # The full path to your galaxy FITS file
psffile = 'gauss25.fits' # Full path to your PSF FITS file
maskfile = 'mask.fits'  # Full path to your mask file
# Do not modify anthing below this comment. These variables are initialised
# here, but are automatically updated by the program. Any changes you
# make here will most likely NOT be used

oldparms = (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
parmnames = []
finalvalues = (0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
# background is recomputed automatically as the median value of the data.

exptime=1.0 # read from the EXPTIME keyword of FITS file, if available
rdnoise=0.0 # read from the RDNOISE keyword of FITS file, if available
gain=1.0    # read from the GAIN keyword of FITS file, if available
ncombine = 1.0 # read from the NCOMBINE keyword of FITS file, if available
vmin=0.0 # minimum value to be shown in progress plots
vmax=0.0 # maximum value to be shown in progress plots

f=fits.open(psffile)
psf=f[0].data
psf=psf/psf.sum()
f.close()

f=fits.open(imagefile)
galaxy = f[0].data # pyfits swaps to python order, second index faster
header = f[0].header
try:
   if (header['exptime']):
     exptime= header['exptime']
   if (header['gain']):
     gain= header['gain']
   if (header['rdnoise']):
     rdnoise= header['rdnoise']
   if (header['ncombine']):
     ncombine= header['ncombine']
 
except:
     print ''
   
f.close()  
nxpts=galaxy.shape[0]
nypts=galaxy.shape[1]
#computedbackground = pylab.median(galaxy.flat)
computedbackground=1000
bulge = np.zeros((nxpts,nypts))
disk = np.zeros((nxpts,nypts))
model_galaxy=np.zeros((nxpts,nypts))
background=computedbackground
if (use_mask):
    f=fits.open(maskfile)
    mask = f[0].data
    f.close()
    # create masked image of the galaxy
    maskedgalaxy = ma.masked_array(galaxy,mask)
    # Set masked values to the median value
    galaxy = ma.filled(maskedgalaxy,fill_value=computedbackground)



