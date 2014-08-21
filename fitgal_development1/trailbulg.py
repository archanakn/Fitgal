
def bfunc(x):
    """ This function gives value of b_n given the Sersic index"""
    return 0.868242*x -0.142058 # Khosroshahi et al. 2000 approximation
    

def trialfunc():
   import numpy as np
   import numpy.random as numra
   import math
   from astropy.io import fits  
   i0b=2000
   re=20
   eb=0.5
   n=4
   i0d=1000
   ed=0.5
   rd=1
   point=0
   bpa=50 * np.pi/180.
   dpa=0 * np.pi/180.
          
   one_minus_eb_sq=(1.0-eb)**2.0
   one_minus_ed_sq=(1.0-ed)**2.0
   nxpts=51
   nypts=51 
   xcenter=26
   ycenter=26 
   xc=xcenter -1.0
   yc=ycenter -1.0 

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
   bulge [xc][yc] = centralpixint/(oversamp*oversamp)
   #if (c.fit_point):
    #    c.bulge [ixc][iyc] = c.bulge [ixc][iyc] + point
   model_galaxy=bulge
   hdu =fits.PrimaryHDU(model_galaxy.astype(np.float32))
   hdu.writeto('test_bulge1.fits',clobber=True)
   print 'completed'
   imag= fits.open('test_bulge.fits') # make 1d array from 2d galaxy array
   galaxy=imag[0].data
   z=galaxy.flat
   numra.seed(120980)# Set random number seed
   zerr = np.sqrt(((abs(numra.poisson(z)-z)*1.0)**2.0+0.0**2.0)) 
   validpix = np.where(zerr != 0.0)
   numpoints = (nxpts*nypts)
   chi2= np.sum(((z[validpix]-model_galaxy.flat[validpix])/zerr[validpix])**2.0)/numpoints 
   print chi2
   return chi2

