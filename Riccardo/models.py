from model import model

const=model()
import math
const._tools.update(__builtins__=None)
md=math.__dict__.copy()
[md.pop(i) for i in math.__dict__ if i.startswith('_')]
const._tools.update(md)
const._tools.update({'int': int})
import numpy
const._tools.update(numpy.__dict__)



const.c='299792458             #m/s #Speed of light'
const.q_e='1.602176462E-19     #C   #Electron charge'
const.m_p='938.271998E+6 *q_e  #kg  #Proton rest mass'
const.m_e='0.510998902E+6*q_e  #kg  #Electron rest mass'
const.h='6.6260693E-34         #Js  #Plank constant'
const.r_p='1.5322E-18          #m   #Proton radius'

lumidata=model()
lumidata.update(const)

lumidata.l_circ    = '26658.883   #m      #LHC circumference'
lumidata.E_p       = '7000  #GeV      #Collision energy'
lumidata.E_coll    = 'E_p*1E9*q_e  #J      #Collision energy'
lumidata.N_b       = '1.15E11     #       #Number of protons per bunch'
lumidata.t_b       = '25E-9       #s      #Bunch spacing'
lumidata.sigma_z   = '11.8E-2     #m      #RMS bunch length'
lumidata.emit_n    = '3.75E-6     #m rad  #Normalized emittance'
lumidata.beta_x    = '55E-2       #m      #Beta x at the IP'
lumidata.beta_y    = '55E-2       #m      #Beta y at the IP'
lumidata.n_b       = '2808        #       #Number of bunches'
lumidata.d_sep     = '9.8         #sigma  #Beam separation in sigma'
lumidata.l_b       = 'c*t_b       #m      #Bunch distance'
lumidata.l_sep     = '130          #m      #Common pipe length'
lumidata.n_par     = 'int(2*l_sep/l_b) #  #Long-range parasitic collisions'
lumidata.f_rev      = 'c/l_circ    #Hz     #Revolution frequency'
lumidata.gamma_r   = 'E_coll/m_p  #            #Relativistic gamma'
lumidata.beta_r    = 'sqrt(1-1/gamma_r**2) #   #Relativistic beta'
lumidata.emit      = 'emit_n/gamma_r/beta_r #m rad   #Emittance'
lumidata.sigma_x   = 'sqrt(emit*beta_x) #m   #RMS beam size at the IP'
lumidata.betastar  = 'beta_x       #m     #Beta at the IP'
lumidata.sigmastar = 'sigma_x       #m     #RMS beam size at the IP'
lumidata.I_b       = 'N_b*n_b*q_e*f_rev  #A      #Beam current'
lumidata.phi_cross   = 'd_sep*sqrt(emit/beta_x) #rad  #Full crossing angle'
lumidata.p_piw     = 'sigma_z*phi_cross/(2*sigma_x) #rad #Piwinski parameter'
lumidata.F_geo     = '1/sqrt(1+p_piw**2)       #   #Geometric factor'
lumidata.L         = 'N_b**2*n_b*f_rev/4/pi/emit/sqrt(beta_x*beta_y)*F_geo #m^2 s^-1 #Luminosity'
lumidata.DeltaQ_head   = 'N_b*r_p/(4*pi*emit*gamma_r) # #Tune shift head on collision'
lumidata.DeltaQ_lr     = '2*n_par*DeltaQ_head/d_sep**2 # #Tune shift long range collision'

lumidata.lstar= '19 #[m] # Distance of the first quadrupole'
lumidata.p_coll    = 'E_coll/c**2*beta_r*c #kg m/s #Momentum at collision'
lumidata.koverg    = 'q_e/p_coll'
lumidata.goverk    = '1/koverg'
lumidata.l1        = '4 #m #length of the first quad'
lumidata.k1        = '1./(lstar*l1)'
lumidata.g1        = 'goverk*k1'
lumidata.bet      = 'beta_x+lstar**2/beta_x'
lumidata.sqbet    = 'sqrt(bet)'
lumidata.sigma     = 'sqbet*sqrt(emit)'
lumidata.aperture     = '33*sigma+7E-3'
lumidata.peakfield = 'g1*aperture/2'

if __name__== '__main__':
  const.graph('const.ps')
  lumidata.graph('lumidata.ps')
  t=model(lumidata)
