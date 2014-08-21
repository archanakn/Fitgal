#!/usr/bin/python
from iminuit import Minuit
from trialbul import *
m=Minuit(main1,i0b=1000.,error_i0b=.1,limit_i0b=(500.,5000.),re=12.,error_re=.1,limit_re=(10.,30.),eb=.3,error_eb=.01,limit_eb=(.1,.8),n=3,error_n=.1,limit_n=(1.,6.),point=0.0,fix_point=True,bpa=30.,error_bpa=.1,limit_bpa=(20.,60.),background=0,fix_background=True)

m.migrad()
