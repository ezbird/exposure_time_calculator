"""
EXPOSURE TIME CALCULATOR
---------------------------
by Khagendra, Priscilla, & Ezra
adapted from Holtz

Purpose:
---------------
  To compute the approximate number of photons you will receive from a given 
  source in a given amount of time for a given observational setup.
 
Details:
---------------
We are solving the "signal" formula, which in practice is difficult because 
we cannot measure all the terms precisely.

S = Tt * integral of (F_lambda/hc/lambda) * a * tel * inst * filt * det * d_lambda

S    is the observed photon flux (the signal in photons/sec)
T    is the telescope collecting area (cm^2)
t    is the integration time, 
a    is the atmospheric transmission
tel  is the efficiency of the telescope
inst is the efficiency of the instrument
filt is the efficiency of the filter
det  is the efficiency of the detector

NOTE: a, and perhaps some of the system efficiencies, are time-variable; 
a is also spatially variable.

Output:
--------------
The final number of photons per second is critical to know in order to estimate  
expected errors and exposure times in observing proposals, observing runs, etc. 
Understanding errors and the uncertainty in the measurement as a function of 
exposure time is absolutely critical.

Classes
--------------
Object
Instrument
Telescope
Mirror
Instrument
Observation

"""
from astropy.modeling import models
from astropy.io import ascii
from astropy import units as u
import astropy.constants as c
from astropy.modeling.models import BlackBody
import numpy as np
import matplotlib.pyplot as plt

class Object:
    '''
    	Astronomical object
        Provide an input spectrum or choose a blackbody.
        Specify temp and magnitude.
        
		Methods: F_nu, F_lambda, photon flux
    '''
    def __init__(self,type='blackbody',temp_eff=10000, mag=0, refmag='V'):
        self.type = type
        self.temp_eff = temp_eff
        self.mag = mag
        self.refmag = refmag
        
    def sed(self,wave):
        # Return SED in F_nu or F_lambda, with units (depending on source of data)
        # NOTE: In the STMAG system, F0, λ = 3.63E - 9ergs/cm2/s/Å, which is the flux of Vega at 5500Å
        if self.type == 'blackbody' :
            the_bb = BlackBody(temperature=self.temp_eff*u.K)
            
            norm= 3.63e-9 * u.erg/u.cm**2 / u.s / u.angstrom / \
                  (the_bb(5500*u.angstrom)*u.sr*c.c/(5500*u.angstrom)**2).to(u.erg/u.cm**2/u.s/u.angstrom)
                  
            return the_bb(wave) * norm*u.sr * 10.**(-0.4*self.mag)
        else :
            raise ValueError('Input type {:s} not yet implemented', self.type)
          
    def photflux(self,wave):
        # Return SED in photon flux
        return (self.f_lambda(wave)*wave/c.h/c.c).to(1/u.cm**2/u.s/u.angstrom)
             
    def photonflux_initial(self) :
        # Returns photon flux from source only
        photflux = (obj.photflux(self.wave)) 
        return photflux
    
    def f_lambda(self,wave):
        # Return SED in F_lambda
        sed=self.sed(wave)
        if sed.unit.is_equivalent(u.erg/u.cm**2/u.s/u.angstrom) :
            # SED is already F_lambda
            return sed
        elif sed.unit.is_equivalent(u.erg/u.cm**2/u.s/u.Hz) :
            # Convert from F_nu to F_lambda
            return (sed*c.c/wave**2).to(u.erg/u.cm**2/u.s/u.angstrom)          
        else: 
            print('SED has incorrect dimensions')
               
    def f_nu(self,wave):
        # Return SED in F_nu
        sed=self.sed(wave)
        if sed.unit.is_equivalent(u.erg/u.cm**2/u.s/u.Hz) :
            # SED is already F_nu
            return sed
        if sed.unit.is_equivalent(u.erg/u.cm**2/u.s/u.angstrom) :
            # convert from Flambda to Fnu
            return (sed*wave**2/c.c).to(1/u.cm**2/u.s/u.Hz)   
        else: 
            print('SED has incorrect dimensions')


class Telescope :
    """ 
    	Telescope
        Given a name (or diameter and mirror array)
        
		Methods: area, throughput
    """
    def __init__(self,name='', diameter=1*u.m, mirrors=['Al']) :
        self.name = name
        if name == 'ARC3.5m' :
            self.diameter=3.5*u.m
            self.mirrors=['Al','Al','Al']
            self.eps=0.4   					# What is this???

        elif name == 'TMO' :
            self.diameter=0.6*u.m
            self.mirrors=['Al','Al']
            self.eps=0.4   # What is this???
        elif name == '' :
            self.diameter=diameter
            self.eps=0.
            self.mirrors=mirrors
        else :
            raise ValueError('Telescope details are unknown.')
    
    def area(self) :
        # Collecting area of the telescope
        return (np.pi*(self.diameter/2)**2*(1-self.eps**2)).to(u.cm**2)
    
    def throughput(self,wave) :
        t=np.ones(len(wave))
        for mir in self.mirrors :
            if type(mir) is float :
                # if mirrors is a float, use it as a constant throughput (bad practice!)
                t *= mir
            else :
                t *= 1
                # Will bring in the mirrors when we have a mirror coating file
                # tmp = Mirror(mir)
                # t *= tmp.reflectivity(wave)
        return t
            
class Mirror :
    """ 
    	Mirror
    	Throughput is a constant percentage at the moment.
    	Default: 0.9 of signal comes through
    	
		TODO: read in file with more detail  a coating name
    """
    def __init__(self,type,const=0.9) :
        self.type = type
        self.const = const
        
    def reflectivity(self,wave) :
        if self.type == 'const' :
            return np.ones(len(wave)) + self.const
            
        else :
            try: 
            	f = open('mirror_coating.txt', 'r')  # We need to re-open the file
            except FileNotFoundError :
                raise ValueError('Coating file not found: {:s}',self.type)
            wav = dat['col1']
            ref = dat['col2']/100.
            interp=scipy.interpolate.interp1d(wav,ref)
            return interp(wave.to(u.nm).value)
    
class Instrument :
    """ 
    	Instrument
    	Throughput is a constant percentage at the moment.
    	Default: 0.8 of signal comes through
    	Also, can read in a filter file which provides 
    	an array of transmission % for each wavelength.
    """
    def __init__(self,name='',efficiency=0.8) :
        self.name = name
        self.efficiency = efficiency
        
        if self.name == 'ARCTIC' :
        	# Imagers ----------------------------------------
        	# We assume the ARCTIC detector for the moment
            self.plate_scale   = 0.228 	 		# arcsec/pixel
            self.binning 	   = 2*2	 	 	# 2x2
            self.readout_noise = 3.7			# Read out noise; 3.7 electrons (slow readout; > 3 seconds)
            self.num_pixels    = 4096*4096   	# 1024 x 1024 pixels
            self.pixel_size    = 15*u.m   		# 15 microns
            self.fov    	   = 7.85   		# arcminutes
        
			# Spectrographs
            if self.name == 'ARCES' :
            	self.readout_noise = 1 # placeholder
        	
    def throughput(self,wave) :
        if self.name == 'ARCTIC' :
            return np.ones(len(wave))*self.efficiency
        else :
            raise ValueError('Do not have this instrument available yet, {:s}',self.name)
        
    def getPixelScale():
    	# Returns the instrument's pixel scale (arcsec/pixel)
    	return self.plate_scale
    	
    def getReadoutNoise():
    	# Returns the instrument's readout noise (electrons)
    	return self.readout_noise
    
    def filter(self,wave,filter='',cent=5500*u.angstrom,width=850*u.angstrom,trans=0.8) :
        if filter=='' :
            out = np.zeros(len(wave))
            out[np.where(( wave>cent-width/2 ) & ( wave<cent+width/2 ))] = trans
            #print(out)
            return out
        else :
            try: 
            	
            	# TODO: issue a warning if the desired wavelengths don't cover the transmission of the filter
            	
            	f = open('APO_filter.txt', 'r')
            	# Loop over lines and extract variables of interest
            	for line in f:
            		line = line.strip()
            		columns = line.split()
            		wavelength = columns[0]
            		trans_percent = columns[1]
            		out = np.zeros(len(wave))
            		print(np.array(trans_percent))
            	return(np.array(trans_percent))
            	
            except FileNotFoundError :
                raise ValueError('No filter file found: {:s} for instrument {:s}',
                                 filter, self.name)


class Atmosphere :
    """ 
    	Earth's atmosphere
    	Throughput is a constant percentage at the moment.
    	Default: 0.8 of signal comes through
    """
    def __init__(self,name='',transmission=0.8) :
        self.name=name
        self.transmission=transmission
    
    def throughput(self,wave) :
        if self.name == '' :
            return np.ones(len(wave))*self.transmission
        else :
            raise ValueError('need to implement atmosphere {:s}',self.name)
            
class Observation :
    """ 
    	The actual observation
    	Send the source's photon flux through the
    	atmosphere, telescope, and instrument.
    	Integrate the photon flux over all wavelengths
    	at the end.
    """
    def __init__(self,obj=None,atmos=None,telescope=None,instrument=None, wave=np.arange(3000,10000,1)*u.angstrom, SNR=None, exp_time=None):
        self.obj 		= obj
        self.atmos 		= atmos
        self.telescope 	= telescope
        self.instrument = instrument
        self.wave		= wave
        self.exp_time	= exp_time
        self.SNR		= SNR
        
    def calc_signal_to_noise(self): # provide exposure time in seconds
        """
        If given exposure time and the object's magnitude, we calculate the SNR
         Just Poisson noise for now: square root of counts
         We need to get the counts/second (N) from the magnitude:
         	m = -2.5 log (N/N0) >>>> 10^(-m/2.5) = N/N0 where N0 is Vega's count rate in the visible
         	obj.mag will eventually be obj.mag = obj.mag + (extinction x airmass)
         Then, SNR = (N x t) / sqrt( N x t + S x p x t + p x R^2 )
         Where S = N x scale^2 x binning^2
        		R = read-out noise (measured in the lab)
        		p = number of pixels in the measuring aperture
        """
        if (self.exp_time == None):
        	raise ValueError('No exposure time provided.')
        
        N 	    	= 1000 * (10**(-obj.mag/2.5)) 		 	# multiply times thousand for vega's base count rate
        plate_scale = inst.plate_scale 	 					# arcsec/pixel
        binning 	= inst.binning 	 	 					# 2x2
        R 			= inst.readout_noise		 			# Read out noise (slow readout; > 3 seconds)
        p  			= inst.num_pixels
        S 	    	= N * plate_scale**2 * binning
        
        print("N: ",N,"plate_scale: ",plate_scale,"binning: ",binning,"R: ",R,"p: ",p,"S: ",S)
        SNR     = (N * self.exp_time) / np.sqrt( (N*self.exp_time) + (S*p*self.exp_time) + (p*R**2) )
        return SNR
    
    def calc_exposure_time(self):
        # If given SNR value and the object's magnitude, we calculate the need exposure time
        # We can solve the SNR equation for t, rearranging to make a quadratic.
        
        if (self.SNR == None):
        	raise ValueError('No signal to noise provided.')
        	
        background_noise = 1 # we'll add background noise later
        N 	    	= 1000 * (10**(-obj.mag/2.5)) 		 	# multiply times thousand for vega's base count rate
        plate_scale = inst.plate_scale 	 					# arcsec/pixel
        binning 	= inst.binning 	 	 					# 2x2
        R 			= inst.readout_noise		 			# Read out noise (slow readout; > 3 seconds)
        p  			= inst.num_pixels
        S 	    	= N * plate_scale**2 * binning
        
        a = N**2
        b = (self.SNR^2) * (N + (p * S) )
        c = (self.SNR^2) * p * R**2

        t = abs( ( -b + np.sqrt(b**2 - 4*a*c) ) / ( 2*a ) )
        
        return t
    
    def photonflux(self) :
        # Returns photon flux
        flux = (obj.photflux(self.wave)* 
                self.atmos.throughput(self.wave) *
                self.telescope.area()*self.telescope.throughput(self.wave)*
                self.instrument.throughput(self.wave)*self.instrument.filter(self.wave) ) 
        return flux
    
    def final_counts_per_second(self) :
        # Returns integral of final photon flux
        return np.trapz(self.photonflux(),self.wave)
        
        
        
# Array of wavelengths to observe
wave  = np.arange(3400,6000,5)*u.angstrom

# Initiative our 4 classes
obj   = Object( type='blackbody', temp_eff=7000, mag=10 ) # Specify an object with temp and magnitude
tel   = Telescope( name='ARC3.5m' )   	  				  # Specify a telescope with size and num of mirrors
inst  = Instrument( name='ARCTIC' )					  	  # Specify the instrument's efficiency
atmos = Atmosphere( transmission=0.8 )					  # Specify the efficiency of the atmosphere

# Call our observation class with our 4 situation classes
# Either a SNR or an exposure time should be provided.
obs   = Observation(obj=obj,atmos=atmos,telescope=tel,instrument=inst,wave=wave,exp_time=100,SNR=1000)

# Throughput of each component, each reducing the final photons per second received
# ---------------------------------------------------------------------------------
# 1.) Photons from the object itself
photflux=obj.photflux(wave)
#plt.plot(wave,photflux, label='Object only', linewidth=4)

# 2.) Photons from the object through the telescope
photflux*=tel.throughput(wave)*tel.area()
plt.plot(wave,photflux, label='Obj + Tel', linewidth=4)

# 3.) Photons from the object through the telescope through the atmosphere
photflux*=atmos.throughput(wave)
plt.plot(wave,photflux, label='Obj + Tel + Atm', linewidth=4)

# 4.) Photons from the object through the telescope through the atmosphere 
#     through the instrument
photflux*=inst.throughput(wave)
plt.plot(wave,photflux, label='Obj + Tel + Atm + Inst', linewidth=4)

# 5.) Photons from the object through the telescope through the atmosphere 
#     through the instrument and its filter(s)
photflux*=inst.filter(wave,trans=0.8,filter='')
plt.plot(wave,photflux, label='Obj + Tel + Atm + Inst + Filter', linewidth=4)

# 6.) With an imaging instrument, we want to integrate over the spectrum 
#     to get the total signal in photons per second
print(obs.final_counts_per_second())

# Print signal to noise and/or exposure time
print("Calculated SNR: ",obs.calc_signal_to_noise(), " from exp time of 100s")
print("Calculated exposure time: ",obs.calc_exposure_time(), " from SNR of 1000")

plt.legend(loc="upper left")
plt.ylabel('Counts')
plt.xlabel('Wavelength (A)')
plt.title('Expected Signal for ARC 3.5m')
plt.show()
