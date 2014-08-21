import numpy as np
import math
from astropy.io import fits
import numpy.random as numra
def main1(i0b,re,eb,n,point,bpa,background):    
    
    parms=np.zeros(11)
    parms[0]=i0b
    parms[1]=re
    parms[2]=eb
    parms[3]=n
    parms[7]=point
    parms[8]=bpa
    parms[10]=background
    generate_model_galaxy(parms)
    imag1=fits.open('test_bulge1.fits')
    model_galaxy=imag1[0].data
    
    imag= fits.open('test_bulge.fits') # make 1d array from 2d galaxy array
    galaxy=imag[0].data
    z= galaxy.flat # make 1d array from 2d galaxy array
    numra.seed(120980)# Set random number seed
    #zerr = np.sqrt(1+0.0*((abs(numra.poisson(z)-z)*c.gain)**2.0+c.rdnoise**2.0)) 
    zerr = np.sqrt((abs(numra.poisson(z)-z)*1.0)**2.0+0.0**2.0) 
   #  print zerr
    
   
    validpix = np.where(zerr != 0.0)
    #if (c.use_mask):
    # numpoints = ma.count(c.maskedgalaxy)
    #else:
    #numpoints = 51 * 51
    
    #chi2= np.sum(((z[validpix]-model_galaxy.flat[validpix])/zerr[validpix])**2.0)/numpoints 
    #print chi2
    avg=np.mean(galaxy-model_galaxy)
    return avg   

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
    point = parms[4]
    bpa = parms[5] * np.pi/180.
    background= parms[6]    
    one_minus_eb_sq=(1.0-eb)**2.0
    nxpts=51
    nypts=51 
    xc=25
    yc=25 
    ixc=int(xc)
    iyc=int(yc)        
    x = np.reshape(np.arange(nxpts*nypts),(nxpts,nypts)) / nypts
    x = x.astype('float32')
    y = np.reshape(np.arange(nxpts*nypts),(nxpts,nypts)) % nypts
    y = y.astype('float32')
        # r is the radius parameter
    xsq=((x-xc)* math.cos(bpa) + (y-yc) * math.sin(bpa))**2.0
    ysq=((xc-x) * math.sin(bpa) + (y-yc) * math.cos(bpa))**2.0
    r=np.sqrt(xsq + ysq/one_minus_eb_sq)
    bulge = i0b * 10.0**(-2.303*(0.868242*n-0.142058)*(r/re)**(1.0/n))
         
        #compute central bulge pixel by oversampling
    centralpixint=0.0
    oversamp=5.0
    for x in np.arange(xc-2.0/oversamp,xc+2.0/oversamp+0.001,1.0/oversamp): 
         for y in np.arange(yc-2.0/oversamp,yc+2.0/oversamp+0.001,1.0/oversamp):
               xsq=((x-xc)* math.cos(bpa) + (y-yc) * math.sin(bpa))**2.0
               ysq=((xc-x) * math.sin(bpa) + (y-yc) * math.cos(bpa))**2.0
                
               r=math.sqrt(xsq + ysq/one_minus_eb_sq )
               try:
                   centralpixint += i0b * 10.0**(-2.303*bfunc(n)*(r/re)**(1.0/n)) 
               except:
                   print i0b,n,bfunc(n),r,re
                   print "Unphysical central pixel parameters, contact YW!"
    bulge [ixc][iyc] = centralpixint/(oversamp*oversamp)
    #if (fit_point):
    #     bulge [ixc][iyc] = c.bulge [ixc][iyc] + point
        


    #if(do_convolve):
    #    c.model_galaxy= conv.convolve2d(c.bulge+c.disk,c.psf,mode='same',boundary='symm')
    #else:
    model_galaxy=bulge
    #if (c.use_mask):
    #    maskedmodel = ma.masked_array(c.model_galaxy,mask=c.mask)
    #    c.model_galaxy =maskedmodel.filled(c.computedbackground)
       
    
    hdu =fits.PrimaryHDU(model_galaxy.astype(np.float32))
    hdu.writeto('test_bulge1.fits',clobber=True)
    #print 'completed'
    #return model_galaxy
