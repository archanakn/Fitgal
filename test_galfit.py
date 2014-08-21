################################################
# Makes up fits images from the bulge and disk 
# derived parameters by GALFIT, i.e. as given in 
# a "fit.log" output file, e.g.
#
################################################
import numpy, pyfits, os, math, scipy.integrate
import pyfits

def get_fit_log_sky():

 f = open('fit.log','r')

 while 1:
   line = f.readline()
   if line == "": break
   word = line.split()   
   if len(word) < 2: continue
   if word[0] == "sky" : 
      skyval = float(word[4])
 f.close()
 return (skyval)

def get_fit_log_size():

 f = open('fit.log','r')
 while 1:
   line = f.readline()
   if line == "": break
   word = line.split()   
   if len(word) < 2: continue
#
# get image size
#
   if word[0] == "Input" : 
      k1 = word[3].find('test1.fits[1:')
      k2 = word[3].find(',1:')
      nx = int(word[3][k1+11:k2])
      k1 = word[3].find(',1:')
      ny = int(word[3][k1+3:-1])
 f.close()
 return (nx,ny)


def get_fit_log_sersic():

 f = open('fit.log','r')

 while 1:
   line = f.readline()
   if line == "": break
   word = line.split()   
   if len(word) < 2: continue
   if word[0] == "sersic" : 
      k1 = word[2].find('(')
      k2 = word[2].find(',')
      ser_xcen = float(word[2][k1+1: k2])
      k1 = 0
      k2 = word[3].find(')')
      ser_ycen = float(word[3][k1: k2])
      ser_mag = float(word[4])
      ser_re  = float(word[5])
      ser_n   = float(word[6])
      ser_boa = float(word[7])
      ser_pa  = float(word[8])
 f.close()

 return (ser_xcen, ser_ycen, ser_mag, ser_re, ser_n, ser_boa, ser_pa)


def get_fit_log_expdisk():

 f = open('fit.log','r')

 while 1:
   line = f.readline()
   if line == "": break
   word = line.split()   
   if len(word) < 2: continue
   if word[0] == "expdisk" : 
      k1 = word[2].find('(')
      k2 = word[2].find(',')
      exp_xcen = float(word[2][k1+1: k2])
      k1 = 0
      k2 = word[3].find(')')
      exp_ycen = float(word[3][k1: k2])
      exp_mag = float(word[4])
      exp_re  = float(word[5])
      exp_boa = float(word[6])
      exp_pa  = float(word[7])
 f.close()

 return (exp_xcen, exp_ycen, exp_mag, exp_re, exp_boa, exp_pa)


def gammln(xx):
#  Logarithm of the gamma function.
        gammln_cof = numpy.array([76.18009173, -86.50532033, 24.01409822,
	                          -1.231739516e0, 0.120858003e-2, -0.536382e-5])
        gammln_stp = 2.50662827465
	x = xx - 1.
	tmp = x + 5.5
	tmp = (x + 0.5)*math.log(tmp) - tmp
	ser = 1.
	for j in range(6):
		x = x + 1.
		ser = ser + gammln_cof[j]/x
	return tmp + math.log(gammln_stp*ser)
#====================================================

def get_b_n(m):
#   find sersic b_n coefficient
#   this is more accurate than the usual b_n = 1.9992*n_sersic - 0.3271
    if m == 0.0: return -0.3271
    b_n = 2.0*m - 1.0/3.0 + 4.0/m/405.0 + 46.0/m/m/25515.0 \
                + 131.0/m/m/m/1148175.0  \
                - 2194697.0/m/m/m/m/30690717750.0
    return b_n

def integrate_pixel_expdisk(x1,x2,y1,y2):
    global yy1, yy2
    yy1 = y1
    yy2 = y2
    return scipy.integrate.dblquad(expdisk,x1,x2, g, h)

def expdisk(y,x):
#  note array order is y, x!!!!
    global xcen, ycen
    global sin_disk, cos_disk, r_disk, b_o_a_disk
    dx = x - xcen
    dy = y - ycen   
    xx = dx*cos_disk + dy*sin_disk
    yy = (-1.0*dx*sin_disk + dy*cos_disk)/b_o_a_disk
    r = math.sqrt(xx*xx + yy*yy)
    return math.exp(-1.0*r/r_disk)

def integrate_pixel_sersic(x1,x2,y1,y2):
    global yy1, yy2
    yy1 = y1
    yy2 = y2
    return scipy.integrate.dblquad(sersic,x1,x2, g, h)

def sersic(y,x):
#  note array order is y, x!!!!
    global xcen, ycen
    global n_sersic, sin_sersic, cos_sersic, r_sersic, b_o_a_sersic
    dx = x - xcen
    dy = y - ycen   
    xx = dx*cos_sersic + dy*sin_sersic
    yy = (-1.0*dx*sin_sersic + dy*cos_sersic)/b_o_a_sersic
    r = math.sqrt(xx*xx + yy*yy)
    return math.exp(-1.*b_n*((r/r_sersic)**(1./n_sersic) - 1.0))

def g(x):
    global yy1
    return yy1

def h(x):
    global yy2
    return yy2

def get_zp(fits):
  base = pyfits.open(fits)
  head = base[0].header
  zeropt = head['MAGZP']
  return zeropt 

#------------------------------------------------------------------------
global xcen, ycen
global sin_disk, cos_disk, r_disk, b_o_a_disk
global n_sersic, sin_sersic, cos_sersic, r_sersic, b_o_a_sersic

hdu=pyfits.open('test1.fits')
raw_head=hdu[0].header
hdu.close()

zeropt = get_zp('test1.fits')
(nx,ny) = get_fit_log_size()
skycts  = get_fit_log_sky()
sky = numpy.zeros((nx,ny), dtype=numpy.float32)
sky.fill(skycts)
#
#  sersic
#
print "Making Sersic bulge ...",
(ser_xcen, ser_ycen, sersic_mag, sersic_re, ser_n, sersic_b_o_a, ser_pa) =get_fit_log_sersic()

xcen = ser_xcen - 1.0
ycen = ser_ycen - 1.0
b_o_a_sersic = sersic_b_o_a
r_sersic = sersic_re
n_sersic = ser_n
pa_sersic = ser_pa - 90.0
sersic_cts = 10**((zeropt - sersic_mag)/2.5)
area = math.pi*r_sersic*r_sersic*b_o_a_sersic
sb_e_cts = 0.5*sersic_cts/area
sb_e_mag = zeropt - 2.5*math.log10(sb_e_cts)
const = 2.5*math.log10(2.*math.pi)
sb = (sersic_mag + const + 5.0*math.log10(r_sersic))
b_n = get_b_n(n_sersic)
n2 = n_sersic*2.
sersic_const = math.exp(gammln(n2))*n_sersic*math.exp(b_n)/(b_n**n2)
i_sersic = sersic_cts/(sersic_const*2.*math.pi*r_sersic*r_sersic*b_o_a_sersic)
sin_sersic = math.sin(pa_sersic*math.pi/180.)
cos_sersic = math.cos(pa_sersic*math.pi/180.)

bulge = numpy.zeros((nx,ny), dtype=numpy.float32)
for iy in range(ny):
  for ix in range(nx):
     dx = float(ix) - xcen
     dy = float(iy) - ycen
     xx = dx*cos_sersic + dy*sin_sersic
     yy = (-1.0*dx*sin_sersic + dy*cos_sersic)/b_o_a_sersic
     r = math.sqrt(xx*xx + yy*yy)
#  note array order is y, x!!!!
     if r < 1.0*r_sersic: 
         x1 = float(ix) - 0.5
         x2 = float(ix) + 0.5
         y1 = float(iy) - 0.5
         y2 = float(iy) + 0.5
         (value,err)=integrate_pixel_sersic(x1,x2,y1,y2)
         bulge[iy,ix] = i_sersic*value
     else:
         bulge[iy,ix] = i_sersic*math.exp(-1.*b_n*((r/r_sersic)**(1./n_sersic) - 1.0))

print "done"
# 
# expdisk
#
print "Making expdisk ...",
(exp_xcen, exp_ycen, exp_mag, exp_re, exp_boa, exp_pa) = get_fit_log_expdisk()

xcen = exp_xcen - 1.0
ycen = exp_ycen - 1.0
disk_mag = exp_mag
r_disk = exp_re
b_o_a_disk = exp_boa
pa_disk =  exp_pa - 90.
disk_cts = 10**((zeropt - disk_mag)/2.5)
i_disk = disk_cts/(2.*math.pi*r_disk*r_disk*b_o_a_disk)
sin_disk = math.sin(pa_disk*math.pi/180)
cos_disk = math.cos(pa_disk*math.pi/180)
disk = numpy.zeros((nx,ny), dtype=numpy.float32)

for iy in range(ny):
  for ix in range(nx):
     dx = float(ix) - xcen
     dy = float(iy) - ycen
     xx = dx*cos_disk + dy*sin_disk
     yy = (-1.0*dx*sin_disk + dy*cos_disk)/b_o_a_disk
     r = math.sqrt(xx*xx + yy*yy)
#  note array order is y, x!!!!
     if r < 1.0*r_disk: 
         x1 = float(ix) - 0.5
         x2 = float(ix) + 0.5
         y1 = float(iy) - 0.5
         y2 = float(iy) + 0.5
         (value,err)=integrate_pixel_expdisk(x1,x2,y1,y2)
         disk[iy,ix] = i_disk*value
     else:     
         disk[iy,ix] = i_disk*math.exp(-1.0*r/r_disk)

print "done"
total = bulge + sky + disk
raw_head.update('MODEL','all','bulge+disk+sky') 
hdu = pyfits.PrimaryHDU(total,raw_head)
hdulist = pyfits.HDUList([hdu])
if os.path.exists('new.fits'): os.unlink('new.fits')
hdulist.writeto('new.fits')

data = bulge
raw_head.update('MODEL','BULGE','bulge only') 
hdu = pyfits.PrimaryHDU(data,raw_head)
hdulist = pyfits.HDUList([hdu])
if os.path.exists('bulge.fits'): os.unlink('bulge.fits')
hdulist.writeto('bulge.fits')

data = bulge + sky
raw_head.update('MODEL','BandS','bulge and sky') 
hdu = pyfits.PrimaryHDU(data,raw_head)
hdulist = pyfits.HDUList([hdu])
if os.path.exists('bulge_sky.fits'): os.unlink('bulge_sky.fits')
hdulist.writeto('bulge_sky.fits')


data = disk
raw_head.update('MODEL','DISK','disk only') 
hdu = pyfits.PrimaryHDU(data,raw_head)
hdulist = pyfits.HDUList([hdu])
if os.path.exists('disk.fits'): os.unlink('disk.fits')
hdulist.writeto('disk.fits')

data = disk + sky
raw_head.update('MODEL','DandS','disk and sky') 
hdu = pyfits.PrimaryHDU(data,raw_head)
hdulist = pyfits.HDUList([hdu])
if os.path.exists('disk_sky.fits'): os.unlink('disk_sky.fits')
hdulist.writeto('disk_sky.fits')


data = disk/(bulge+disk)
raw_head.update('MODEL','DoverT','Disk over Total counts') 
hdu = pyfits.PrimaryHDU(data,raw_head)
hdulist = pyfits.HDUList([hdu])
if os.path.exists('div.fits'): os.unlink('div.fits')
hdulist.writeto('div.fits')
