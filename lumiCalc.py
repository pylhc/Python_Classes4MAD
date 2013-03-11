from scipy.interpolate import interp1d
from scipy.integrate import quad, dblquad
from pylab import *
from numpy import *
from scipy import stats
from math import *


def getlumi0(bint1,bint2,epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,frev,nbunch,gamma,betarel):
  epsx1 = epsx1n/(gamma*betarel)
  epsy1 = epsy1n/(gamma*betarel)
  epsx2 = epsx2n/(gamma*betarel)
  epsy2 = epsy2n/(gamma*betarel)  
  lumi = bint1*bint2*frev*nbunch/(2*pi*sqrt(epsx1*betax1+epsx2*betax2)*sqrt(epsy1*betay1+epsy2*betay2))*1e-4
  return lumi
  
def getGeomFact(eps1n,beta1,eps2n,beta2,gamma,betarel,sigs1,sigs2,phi):  
  eps1 = eps1n/(gamma*betarel)
  eps2 = eps2n/(gamma*betarel)
  geom = sqrt(1+(sigs1*sigs1+sigs2*sigs2)/(eps1*beta1+eps2*beta2)*tan(phi/2)**2)
  return geom
  
def getBBparam(bint,epsxn,epsyn,betax,betay,gamma,betarel,r0,geom):
  epsx = epsxn/(gamma*betarel)
  epsy = epsyn/(gamma*betarel)
  sigx = sqrt(betax*epsx)
  sigy = sqrt(betay*epsy)
  xi = []
  xi.append(bint*r0*betax/(2*pi*gamma*sigx*geom*(sigx*geom+sigy)))
  xi.append(bint*r0*betay/(2*pi*gamma*sigy*(sigx+sigy*geom)))
  return xi[0]+xi[1]
  
def getHGfact(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2):
  epsx1 = epsx1n/(gamma*betarel)
  epsy1 = epsy1n/(gamma*betarel)
  epsx2 = epsx2n/(gamma*betarel)
  epsy2 = epsy2n/(gamma*betarel)  
  sigx1 = sqrt(epsx1*betax1)
  sigx2 = sqrt(epsx2*betax2) 
  sigy1 = sqrt(epsy1*betay1)
  sigy2 = sqrt(epsy2*betay2)
  tx = 2*(sigx1**2+sigx2**2)/((sigs1**2+sigs2**2)*(sigx1**2/betax1**2+sigx2**2/betax2**2))
  ty = 2*(sigy1**2+sigy2**2)/((sigs1**2+sigs2**2)*(sigy1**2/betay1**2+sigy2**2/betay2**2))
  hgfact = quad(lambda t: 1/sqrt(pi)*exp(-t**2)/sqrt((1+t**2/tx)*(1+t**2/ty)), -inf, inf)
  return hgfact[0]
  
def getReducHGX(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2,phi):
  epsx1 = epsx1n/(gamma*betarel)
  epsy1 = epsy1n/(gamma*betarel)
  epsx2 = epsx2n/(gamma*betarel)
  epsy2 = epsy2n/(gamma*betarel)  
  sigx1 = sqrt(epsx1*betax1)
  sigx2 = sqrt(epsx2*betax2) 
  sigy1 = sqrt(epsy1*betay1)
  sigy2 = sqrt(epsy2*betay2)
  const = sqrt(2.0)*cos(phi/2)**2/(sqrt(pi)*sqrt(sigs1**2+sigs2**2))
  a = 2*sin(phi/2)**2/(sigx1**2+sigx2**2)
  b = 2*cos(phi/2)**2/(sigs1**2+sigs2**2)
  reduc = quad(lambda s: const*exp(-(a/(sqrt(1+s**2/betax1**2)*sqrt(1+s**2/betax2**2))+b)*s**2)/(sqrt(1+s**2/betax1**2)*sqrt(1+s**2/betax2**2)), -inf, inf)
  return reduc[0]
  
def getLumi(lumi0,reduc):
  return lumi0*reduc
  
  
def getIntensity(t,levLumi,xsec,initInt,nip,nbunch):
  rate = nip*levLumi*xsec*1.0e-27/nbunch
  return initInt-rate*t
  
def getPhi(bint1,bint2,epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,frev,nbunch,gamma,betarel,sigs1,sigs2,phi,levLumi):
  phi0 = phi
  lumi = getlumi0(bint1,bint2,epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,frev,nbunch,gamma,betarel)
  reduc = getReducHGX(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2,phi0)
  lumi0 = getLumi(lumi,reduc)
  while(lumi0>levLumi and phi0>0.0):
    phi0 = phi0+1.0e-6
    lumi = getlumi0(bint1,bint2,epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,frev,nbunch,gamma,betarel)
    reduc = getReducHGX(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2,phi0)
    lumi0 = getLumi(lumi,reduc)  
  while(lumi0<levLumi and phi0>0.0):
    phi0 = phi0-1.0e-6
    lumi = getlumi0(bint1,bint2,epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,frev,nbunch,gamma,betarel)
    reduc = getReducHGX(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2,phi0)
    lumi0 = getLumi(lumi,reduc)   
  return phi0
  
  
frev = 11245
nbunch = 2808
nip = 2
r0 = 1.535e-18

energy = 7000e9
pmass = 0.93827231e9
gamma = energy/pmass+1.0
betarel = sqrt(1.0-1.0/(gamma*gamma))

bint1 = 2.2e11
bint2 = 2.2e11

epsx1n = 2.5e-6
epsx2n = 2.5e-6
epsy1n = 2.5e-6
epsy2n = 2.5e-6

betax1 = 0.15
betax2 = 0.15
betay1 = 0.15
betay2 = 0.15

sigs1 = 0.075
sigs2 = 0.075

phi = 590e-6

t = 15*3600
dt = 600

levLumi = 5.0e34
xsec = 100

hgfacti = getHGfact(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2)

fout = open('level.out','w')

phii = phi
lumii = levLumi
bint1i = bint1
bint2i = bint2

for i in range(0,t,dt):
  bint1tmp = getIntensity(dt,lumii,xsec,bint1i,nip,nbunch)
  bint2tmp = getIntensity(dt,lumii,xsec,bint2i,nip,nbunch)
  bint1i = bint1tmp
  bint2i = bint2tmp
  lumi0i = getlumi0(bint1i,bint2i,epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,frev,nbunch,gamma,betarel)
  if(phii>0.0):
    phitmp = getPhi(bint1i,bint2i,epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,frev,nbunch,gamma,betarel,sigs1,sigs2,phii,levLumi)
  else:
    phii = 0.0
  if(phitmp>0.0):
    phii = phitmp
  else:
    phii = 0.0
  geomi = getGeomFact(epsx1n,betax1,epsx2n,betax2,gamma,betarel,sigs1,sigs2,phii)
  hgfacti = getHGfact(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2)
  reduci = getReducHGX(epsx1n,betax1,epsy1n,betay1,epsx2n,betax2,epsy2n,betay2,sigs1,sigs2,phii)
  lumii = getLumi(lumi0i,reduci)
  xi = getBBparam(bint1i,epsx1n,epsy1n,betax1,betay1,gamma,betarel,r0,geomi)
  fout.write(str(i)+"\t"+str(bint1i)+"\t"+str(bint2i)+"\t"+str(lumi0i)+"\t"+str(phii)+"\t"+str(1/geomi)+"\t"+str(hgfacti)+"\t"+str(reduci)+"\t"+str(lumii)+"\t"+str(xi)+"\n")
  
fout.close()
sys.exit(0)