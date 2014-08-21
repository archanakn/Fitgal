#!/usr/bin/env python

import pyfits,os,time,sys,math,glob
import config as c

# parse command line variables
if (len(sys.argv) != 8):
    print "Expected 7 command line arguments. Found "+str(len(sys.argv)-1)+" Using default values for command line arguments\n"
else:
    c.imagefile=sys.argv[1]
    c.psffile = sys.argv[2]
    c.xcenter=  float(sys.argv[3])
    c.ycenter=  float(sys.argv[4])
    if (int(sys.argv[5]) == 1): c.fit_bulge=True
    if (int(sys.argv[6]) == 1): c.fit_disk=True
    if (int(sys.argv[7]) == 1): c.fit_point=True

# Is graphics plotting module available?
try:
    import pylab
    c.have_pylab=True
except:
    print "Plotting library not found"

# import required modules
from numarray import * # numarray provides fast array processing
import random, numarray.random_array,numarray.ma
import numarray.convolve as conv # convolution functions
from hippo import FunctionBase # hippo provides Python-Minuit interface

from fitgalfunc import * # fitgalfunc is a collection of useful functions
      
# Initialise variables
tstart = time.clock() # Start clock to measure time taken for fitting 
c.niter=0 # initialise number of iterations

# Read in the galaxy data and header parameters using pyfits module
f=pyfits.open(c.imagefile)
c.galaxy = f[0].data # pyfits swaps to python order, second index faster
header = f[0].header
if (header.has_key('exptime')):
    c.exptime= header['exptime']
if (header.has_key('gain')):
    c.gain= header['gain']
if (header.has_key('rdnoise')):
    c.rdnoise= header['rdnoise']
if (header.has_key('ncombine')):
    c.ncombine= header['ncombine']
f.close()         
c.nxpts = c.galaxy.shape[0]
c.nypts = c.galaxy.shape[1]

c.computedbackground = pylab.median(c.galaxy.flat)

if (c.use_mask):
    f=pyfits.open(c.maskfile)
    c.mask = f[0].data
    f.close()
    # create masked image of the galaxy
    c.maskedgalaxy = numarray.ma.masked_array(c.galaxy,c.mask)
    # Set masked values to the median value
    c.galaxy = numarray.ma.filled(c.maskedgalaxy,value=c.computedbackground)

# Initialize
c.bulge = numarray.zeros((c.nxpts,c.nypts))
c.disk = numarray.zeros((c.nxpts,c.nypts))

z= c.galaxy.flat # make 1d array from 2d galaxy array

# read in the psf image
f=pyfits.open(c.psffile)
c.psf = f[0].data
f.close()        

# normalise the psf
c.psf = c.psf/c.psf.sum()

# Initialise the arrays that will form the N Tuple for minimization
x = numarray.zeros(c.nxpts*c.nypts)
y = numarray.arange(c.nxpts*c.nypts)
xerr = numarray.zeros(c.nxpts*c.nypts)
yerr = numarray.zeros(c.nxpts*c.nypts)

numarray.random_array.seed(11198205, 120980)# Set random number seed

# compute errors assuming poisson fluctuations in electrons
zerr = numarray.sqrt((abs(numarray.random_array.poisson(z)-z)*c.gain)**2.0+c.rdnoise**2.0)

# Populate x and y arrays
k = 0
for i in range(c.nxpts) :
    for j in range(c.nypts) :
        y[k] = j
        x[k] = i
        k += 1

# Generate the N-Tuple for the input data. Hippodraw needs it in this form

from hippo import NTuple
nt = NTuple () # empty one
nt.addColumn ( 'x', x )
nt.addColumn ( 'y', y )
nt.addColumn ( 'z', z )
nt.addColumn ( 'xerr', xerr )
nt.addColumn ( 'yerr', yerr )
nt.addColumn ( 'zerr', zerr )

# Input data has been read, now fitter function need to configured
from hippo import FitterFactory
factory = FitterFactory.instance()
fitters = factory.names()
fittertouse = fitters[1]
# fitters[1] -> Migrad Chi^2 minimisation
migrad = factory.create (fittertouse)

f = GalaxyFit()

fcn = migrad.getFCN ()
fcn.setFunction ( f ) # f in the function to minimize

# set limits on each parameter setLimits('parameter_name',minval,maxval)
if (fittertouse == 'ChiSq: Minuit(Migrad)'):

    migrad.setLimits('i0b',c.i0bvals[1],c.i0bvals[2])
    migrad.setLimits('re',c.revals[1],c.revals[2])
    migrad.setLimits('eb',c.ebvals[1],c.ebvals[2])
    migrad.setLimits('n',c.nvals[1],c.nvals[2])
    migrad.setLimits('i0d',c.i0dvals[1],c.i0dvals[2])
    migrad.setLimits('rd',c.rdvals[1],c.rdvals[2])
    migrad.setLimits('ed',c.edvals[1],c.edvals[2])
    migrad.setLimits('point',c.pointvals[1],c.pointvals[2])
    migrad.setLimits('bpa',c.bpavals[1],c.bpavals[2])
    migrad.setLimits('dpa',c.dpavals[1],c.dpavals[2])
    migrad.setLimits('background',c.backgroundvals[1],c.backgroundvals[2])
    
    # Set step sizes for each variable
    migrad.setStepSize('i0b',c.i0bvals[3])
    migrad.setStepSize('re',c.revals[3])
    migrad.setStepSize('eb',c.ebvals[3])
    migrad.setStepSize('n',c.nvals[3])
    migrad.setStepSize('i0d',c.i0dvals[3])
    migrad.setStepSize('rd',c.rdvals[3])
    migrad.setStepSize('ed',c.edvals[3])
    migrad.setStepSize('point',c.pointvals[3])
    migrad.setStepSize('bpa',c.bpavals[3])
    migrad.setStepSize('dpa',c.dpavals[3])
    migrad.setStepSize('background',c.backgroundvals[3])

    # set which parameters are fixed
    migrad.setFixedFlags([c.i0bvals[4],c.revals[4],c.ebvals[4],c.nvals[4],c.i0dvals[4],c.rdvals[4],c.edvals[4],c.pointvals[4],c.bpavals[4],c.dpavals[4],c.backgroundvals[4]])

fcn.setDataSource (nt,2,[ 0, 1, 2, -1, -1, 5]) # tell fn what data to use

fcn.setUseErrors(True) # Use Errors (True) or not (False)

result=migrad.minimize() # Perform the minimization using Minuit

# Minimization should now be over, start to generate outputs

c.finalvalues=f.getParameters() # read final fit parameters into list

# remove old files if they exist
for myfile in ['model.fits','bulge.fits','disk.fits','residual.fits']:
    if os.access(myfile,os.F_OK):
        os.remove(myfile)

# Write Model galaxy image
hdu = pyfits.PrimaryHDU(c.model_galaxy.astype(Float32))
hdu.writeto('model.fits')

# Write bulge
c.bulge= conv.convolve2d(c.bulge,c.psf,fft=0,mode='constant')
hdu = pyfits.PrimaryHDU(c.bulge.astype(Float32))
hdu.writeto('bulge.fits')

# Write disk
c.disk = conv.convolve2d(c.disk,c.psf,fft=0,mode='constant')
hdu = pyfits.PrimaryHDU(c.disk.astype(Float32))
hdu.writeto('disk.fits')

# Write residual image defined as in Galfit
residual = c.galaxy-c.model_galaxy
hdu = pyfits.PrimaryHDU(residual.astype(Float32))
hdu.writeto('residual.fits')

# Write fit results as a CSV format file
resultfile = 'fitgal.csv'
outfile = open(resultfile,'w')
outfile.write(str(f.parmNames())[1:-1]+'\n')
outfile.write(str(f.getParameters())[1:-1]+'\n')
outfile.close()

# Generate plot and write as PNG and PS files
if (c.have_pylab and c.make_plots):
    write_html_norefresh()
    pylab.subplot(223)
    pylab.savefig('fitgal.ps')
    # Generate animated GIF for the fitting process
    if (c.make_movie):
        os.system("convert -compress LZW -adjoin -delay 400 iter????.png fitgal-animation.gif")
        files=glob.glob('iter????.png')
        for file in files:
            os.remove(file)

# Compute time taken and chi^2
tend=time.clock() # CPU time
timetaken= tend - tstart # tstart was defined earlier
validpix = where(zerr != 0.0)

if (c.use_mask):
    numpoints = numarray.ma.count(c.maskedgalaxy)
else:
    numpoints = (c.nxpts*c.nypts)
    
chi2 = numarray.sum(((z[validpix]-c.model_galaxy.flat[validpix])/zerr[validpix])**2.0)/numpoints

print "Migrad fitting completed in "+str(c.niter)+" iterations with reduced chi^2 "+str(chi2)+" in "+str(timetaken)+" seconds"
