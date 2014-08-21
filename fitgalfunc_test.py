#!/usr/bin/env python
import sys
import random,numpy.random,math
import scipy.signal as conv
import config_test as c
import numpy as np
from iminuit import Minuit 
import numpy.ma as ma
from time import sleep
from astropy.io import fits
try: # Check if plotting routines are available
    import pylab
    c.have_pylab=True
    pylab.ion()
except:
    print "plotting library not found, plots will NOT be generated!"

# Transform to Python conventions (0 indexed, second index varies faster)
c.xc = c.ycenter-1.0
c.yc= c.xcenter - 1.0
def valueAt(i0b,re,eb,n,i0d,rd,ed,point,bpa,dpa,background):
       
       c.parmnames = ( [ 'i0b', 're', 'eb','n','i0d','rd','ed','point','bpa','dpa','background'] );
       ix=c.xcenter
       iy=c.ycenter
       parms=np.zeros(11)
       parms[0]=i0b
       parms[1]=re
       parms[2]=eb
       parms[3]=n
       parms[4]=i0d
       parms[5]=rd
       parms[6]=ed
       parms[7]=point
       parms[8]=bpa
       parms[9]=dpa
       parms[10]=background
       #print 'parms'
       #print parms
       #print 'oldparms'
       #print c.oldparms
       #sys.stdout.flush()
       #sleep(0.5) 
       if (np.any(parms != c.oldparms)):
           generate_model_galaxy(parms)
           str1=str(i0b)+''+str(re)+' '+str(eb)+' '+str(n)+' '+str(i0d)+' '+str(rd)+' '+str(ed)+' '+str(point)+' '+str(bpa)+' '+str(dpa)+' '+str(background) 
           
           
           
           c.niter += 1
           c.oldparms = parms
           if (c.niter % c.plotiter == 1):
               print "Iteration:",c.niter,parms
               pylab.ion()  
               plot_fitstatus(str1)
               if (c.have_pylab and c.make_plots):
                   write_html()
       #sys.stdout.flush()
       #sleep(0.5)
       return c.model_galaxy[ix][iy]        


def bfunc(x):
    """ This function gives value of b_n given the Sersic index"""
    return 0.868242*x -0.142058 # Khosroshahi et al. 2000 approximation
    
def generate_model_galaxy(parms):
    
    """Function to generate the model galaxy with various components"""
    make_bulge=False
    make_disk=False
    import math
    i0b=parms[0]
    re=parms[1]
    eb=parms[2]
    n=parms[3]
    i0d = parms[4]
    rd = parms[5]
    ed = parms[6]
    point = parms[7]
    bpa = parms[8] * np.pi/180.
    #print 'bpa' 
    #print bpa 
    dpa =  parms[9] * np.pi/180.
    #print 'dpa'
    #print dpa
    oldbpa = c.oldparms[8] *np.pi/180. # convert to radians
    olddpa = c.oldparms[9] *np.pi/180. # convert to radians
    background= parms[10]
    c.background = parms[10]
    #print 'background'
    #print background
    #print 'c.background'
    #print c.background

    if (bpa!=oldbpa or i0b!=c.oldparms[0] or re!=c.oldparms[1] or eb != c.oldparms[2] or n!=c.oldparms[3] or point!=c.oldparms[7]):
        make_bulge=True
    
    
    
    if (dpa!=olddpa or i0d!=c.oldparms[4] or rd!=c.oldparms[5] or ed != c.oldparms[6]):
        make_disk=True
        
        
    one_minus_eb_sq=(1.0-eb)**2.0
    one_minus_ed_sq=(1.0-ed)**2.0
    
    ixc=int(c.xc)
    iyc=int(c.yc)
    
   
    if (make_bulge and c.fit_bulge):
        
        x = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) / c.nypts
        x = x.astype('float32')
        y = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) % c.nypts
        y = y.astype('float32')
        # r is the radius parameter
       # math.cos(dpa)=math.cos(bpa)
        #math.sin(dpa)=math.sin(bpa)
        xsq=((x-c.xc)* math.cos(bpa) + (y-c.yc) * math.sin(bpa))**2.0
        ysq=((c.xc-x) * math.sin(bpa) + (y-c.yc) * math.cos(bpa))**2.0
        r=np.sqrt(xsq + ysq/one_minus_eb_sq)
        c.bulge = i0b * 10.0**(-2.303*(0.868242*n-0.142058)*(r/re)**(1.0/n))
        
        #compute central bulge pixel by oversampling
        centralpixint=0.0
        oversamp=5.0
        for x in np.arange(c.xc-2.0/oversamp,c.xc+2.0/oversamp+0.001,1.0/oversamp): 
            for y in np.arange(c.yc-2.0/oversamp,c.yc+2.0/oversamp+0.001,1.0/oversamp):
                xsq=((x-c.xc)* math.cos(bpa) + (y-c.yc) * math.sin(bpa))**2.0
                ysq=((c.xc-x) * math.sin(bpa) + (y-c.yc) * math.cos(bpa))**2.0
                
                r=math.sqrt(xsq + ysq/one_minus_eb_sq )
                try:
                    centralpixint += i0b * 10.0**(-2.303*bfunc(n)*(r/re)**(1.0/n)) 
                except:
                    print i0b,n,bfunc(n),r,re
                    print "Unphysical central pixel parameters, contact YW!"
        c.bulge [ixc][iyc] = centralpixint/(oversamp*oversamp)
        if (c.fit_point):
            c.bulge [ixc][iyc] = c.bulge [ixc][iyc] + point
        

    if (make_disk and c.fit_disk):
        x = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) / c.nypts
        x = x.astype('float32')
        y = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) % c.nypts
        y = y.astype('float32')
        #math.cos(dpa)=math.cos(dpa)
        #math.sin(dpa)=math.sin(dpa)
        xsq=((x-c.xc)* math.cos(dpa) + (y-c.yc) *math.sin(dpa))**2.0
        ysq=((c.xc-x) * math.sin(dpa) + (y-c.yc) *math.cos(dpa))**2.0
        r=np.sqrt(xsq + ysq/one_minus_ed_sq)
        c.disk  = i0d * np.exp(-1.0 * np.array(r/rd))
    if(c.do_convolve):
        c.model_galaxy= conv.convolve2d(c.bulge+c.disk+c.background,c.psf,mode='same',boundary='symm')
    if (c.use_mask):
        maskedmodel = ma.masked_array(c.model_galaxy,mask=c.mask)
        c.model_galaxy =maskedmodel.filled(c.computedbackground)
        #print maskedmodel
        #print c.model_galaxy

    

def plot_fitstatus(str1):
   

    """Function to plot the progress of the fitting process"""
    #if (c.niter == 1):
    pylab.ion()
    pylab.clf()
    pylab.subplot(221)
    image1=pylab.imshow(np.rot90(np.swapaxes(c.galaxy,0,1)))
    image1.autoscale() # This line is necessary to set image1.norm.v*
    c.vmin= image1.norm.vmin
    c.vmax= image1.norm.vmax
    pylab.title(c.imagefile)
    
        
    pylab.subplot(222)
    pylab.ion()
    pylab.imshow(np.rot90(np.swapaxes(c.model_galaxy,0,1)),aspect='equal',vmin=c.vmin,vmax=c.vmax)
    if (c.niter == 1):
        pylab.colorbar()
    pylab.title('Model Galaxy, Iteration='+str(c.niter))
    pylab.subplot(223)
    pylab.ion()
    residual = np.zeros((c.nxpts,c.nypts))
    residual2 = np.zeros((c.nxpts,c.nypts))
    validpix = np.where ((c.galaxy-c.background) !=0.0)
    #print validpix
    residual[validpix] = c.galaxy[validpix]-c.model_galaxy[validpix]
    residual2[validpix] = c.galaxy[validpix]-c.model_galaxy[validpix]
    notvalidpix =  np.where ((c.galaxy-c.background)==0.0)
    residual2[notvalidpix] = 0.0
    histrange=(c.vmax-c.vmin)/5.0
    pylab.imshow(np.rot90(np.swapaxes(residual2,0,1)),aspect='equal',vmin=-1.0*histrange,vmax=histrange)
    if (c.niter == 1):
        pylab.colorbar()
    pylab.title('Residual = Real-Model')
    
    #pylab.subplot(224)
    # the histogram of the data
    pylab.text(70,25,str1)
    #pylab.hist(residual,bins=np.arange(-1.0*histrange,histrange),hold=False)
    #pylab.title('Difference Histogram')
    if (c.make_movie):
        iternum =  "%04d" % (c.niter)
        pylab.savefig('iter'+iternum+'.png')
        
       #  sys.stdout.flush()
       # sleep(0.10)
    
def write_html():
    #plot_fitstatus('hi')
    #pylab.savefig('fitgal.png')
    outfile = open('fitgal.html','w')
    outfile.write('<HTML><HEAD>')
    outfile.write('<META HTTP-EQUIV="Refresh" CONTENT="10; URL='+c.webdir+'fitgal.html">')
    outfile.write('</HEAD><BODY>')
    outfile.write('<CENTER><IMG SRC="fitgal.png"></CENTER>')
    outfile.write('Iteration:'+str(c.niter)+'<BR>')
    outfile.write('<B>'+str(c.parmnames)+'<BR>'+str(c.oldparms)+'</B>')
    outfile.write('</BODY></HTML>')
    outfile.close()


def write_html_norefresh():
    #plot_fitstatus('hi')
    #pylab.savefig('fitgal.png')
    f = open('fitgal.html','w')
    f.write('<HTML><HEAD>')
    f.write('</HEAD><BODY>')
    f.write('<CENTER><H1> All done!</H1><IMG SRC="fitgal.png"></CENTER>')
    f.write('<HR><H1>Fitting process completed </H1>')
    # FIX D/B is correct only for n=4, correct this
    c.dbyb = math.exp(7.67) * (c.finalvalues[4]/c.finalvalues[0]) * (c.finalvalues[5]/c.finalvalues[1])**2.0
    if (c.fit_bulge and c.fit_disk):
	f.write('<H1> D/B ratio:'+str(c.dbyb)+'<BR>')
    if (c.fit_bulge):
        f.write('<H2>Bulge</H2>')
        f.write('Central intensity (DN): '+str(c.finalvalues[0])+'<BR>')
        f.write('Half light radius (pixels): '+str(c.finalvalues[1])+'<BR>')
        f.write('Sersic Index n: '+str(c.finalvalues[3])+'<BR>')
        f.write('Ellipticity: '+str(c.finalvalues[2])+'<BR>')
    if (c.fit_disk):
        f.write('<H2>Disk</H2>')
        f.write('Central Intensity (DN): '+str(c.finalvalues[4])+'<BR>')
        f.write('Scale length (pixels): '+str(c.finalvalues[5])+'<BR>')
        f.write('Ellipticity: '+str(c.finalvalues[6])+'<BR>')
    if (c.fit_point):
        f.write('<H2>Point</H2>')
        f.write('Point intensity (DN): '+str(c.finalvalues[7])+'<BR>')
    f.write('<H2>Other parameters</H2>')
    f.write('Bulge Position angle (degrees): '+str(c.finalvalues[8])+'<BR>')
    f.write('Disk Position angle (degrees): '+str(c.finalvalues[9])+'<BR>')
    f.write('Background (DN): '+str(c.finalvalues[10])+'<BR>')
    f.write('<H1>Download output files</H1>')
    f.write('<UL><LI>Best fit <A HREF="model.fits">model</A> galaxy image (FITS)')
    f.write('<LI>Best fit <A HREF="bulge.fits">bulge</A> component image (FITS)')
    f.write('<LI>Best fit <A HREF="disk.fits">disk</A> component image (FITS)')
    f.write('<LI>Best fit <A HREF="residual.fits">residual</A>  image (FITS)')
    f.write('<LI>Best fit parameters <A HREF="fitgal.csv">list</A> (CSV)')
    f.write('<LI>Graphical fit output summary (<A HREF="fitgal.png">PNG</A>,<A HREF="fitgal.ps">EPS</A>)')
    
    f.write('</UL></BODY></HTML>')
    f.close()         

        


def make_test_galaxy():
     i0b=2000
     re=4
     eb=0.2
     n=4
     i0d=1000
     ed=0.5
     rd=4
     point=0
     bpa=30 * np.pi/180.
     dpa=30 * np.pi/180.
     background=1300        
          
     one_minus_eb_sq=(1.0-eb)**2.0
     one_minus_ed_sq=(1.0-ed)**2.0
    
     ixc=int(c.xc)
     iyc=int(c.yc)

     x = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) / c.nypts
     x = x.astype('float32')
     y = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) % c.nypts
     y = y.astype('float32')
     # r is the radius parameter
     xsq=((x-c.xc)* math.cos(bpa) + (y-c.yc) * math.sin(bpa))**2.0
     ysq=((c.xc-x) * math.sin(bpa) + (y-c.yc) * math.cos(bpa))**2.0
     r=np.sqrt(xsq + ysq/one_minus_eb_sq)
     c.bulge = i0b * 10.0**(-2.303*(0.868242*n-0.142058)*(r/re)**(1.0/n))
        
     #compute central bulge pixel by oversampling
     centralpixint=0.0
     oversamp=5.0
     for x in np.arange(c.xc-2.0/oversamp,c.xc+2.0/oversamp+0.001,1.0/oversamp): 
        for y in np.arange(c.yc-2.0/oversamp,c.yc+2.0/oversamp+0.001,1.0/oversamp):
              xsq=((x-c.xc)* math.cos(bpa) + (y-c.yc) * math.sin(bpa))**2.0
              ysq=((c.xc-x) * math.sin(bpa) + (y-c.yc) * math.cos(bpa))**2.0
                
              r=math.sqrt(xsq + ysq/one_minus_eb_sq )
              try:
                  centralpixint += i0b * 10.0**(-2.303*bfunc(n)*(r/re)**(1.0/n)) 
              except:
                   print i0b,n,bfunc(n),r,re
                   print "Unphysical central pixel parameters, contact YW!"
     c.bulge [ixc][iyc] = centralpixint/(oversamp*oversamp)
     
     c.bulge [ixc][iyc] = c.bulge [ixc][iyc] + point

     x = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) / c.nypts
     x = x.astype('float32')
     y = np.reshape(np.arange(c.nxpts*c.nypts),(c.nxpts,c.nypts)) % c.nypts
     y = y.astype('float32')
     xsq=((x-c.xc)* math.cos(dpa) + (y-c.yc) *math.sin(dpa))**2.0
     ysq=((c.xc-x) * math.sin(dpa) + (y-c.yc) *math.cos(dpa))**2.0
     r=np.sqrt(xsq + ysq/one_minus_ed_sq)
     c.disk  = i0d * np.exp(-1.0 * np.array(r/rd))
     c.model_galaxy= conv.convolve2d(c.bulge+c.disk+c.background,c.psf,mode='same',boundary='symm')
     if (c.use_mask):
        maskedmodel = ma.masked_array(c.model_galaxy,mask=c.mask)
        c.model_galaxy =maskedmodel.filled(c.computedbackground) 

     hdu =fits.PrimaryHDU(c.model_galaxy.astype(np.float32))
     hdu.writeto('model_galaxy.fits')
     print 'completed'
