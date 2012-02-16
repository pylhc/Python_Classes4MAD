

import sys, gzip
try:
    from numpy import *
except ImportError:
    try:
        from Numeric import *
        from LinearAlgebra import *
    except ImportError: sys.exit()

from string import split
from string import replace

def replac(nm): #-- added aug12_08 calaga
    nm=replace(nm,"-","_")
    nm=replace(nm,".","_")
    return nm


#########################
class YASPtwiss:
#########################
    "Read orbit data from YASP output in the same way as metaclass"
    def __init__(self, filename=''): 


        self.HNAME=[]
        self.HCO=[]
        self.HCOerr=[]
        self.Hindx={}
        self.VNAME=[]
        self.VCO=[]
        self.VCOerr=[]
        self.Vindx={}
        self.HcorrNAME=[]
        self.Hcorr=[]
        self.VcorrNAME=[]
        self.Vcorr=[]

        if filename=='':
            return

        
        if '.gz' in filename:
            f=gzip.open(filename, 'rb')
        else:
            f=open(filename, 'r')



        #fspl=filename.split('/')[-1].split('.')
        #self.exciter=fspl[1]+fspl[2]
        #self.plane=fspl[3]
        #if "minus" in fspl[4]:
        #  self.str=-1.0
        #else:
        #  self.str=1.0
        


        corr=0
        count=0
        for line in f:
            cl=line.split()
            if ("@ " in line) :
                label=replac(cl[1])
		try:
                	exec "self."+label+"= '"+str((cl[3]))+"'"
		except:
			print "Problem parsing:", line
			#print "Ignored"

            if "# CORRECTOR" in line:
                corr=1

            #if ("* " in line) :
            #    alllabels=split(line)
            #    print "alllabels",len(alllabels)
            #    for j in range(1,len(alllabels)):
            #        exec "self."+alllabels[j]+"= []"



            
            if corr==0 and cl[1]=='H' and cl[8]=="OK":
                
                self.HCO.append(float(cl[3])*1e-3)
                self.HCOerr.append(0.0)
                #self.Hindx[replace(cl[0], '"', '')]=len(self.HNAME)-1
                #self.Hindx[replace(cl[0], '"', '').upper()]=len(self.HNAME)-1
                
                self.HNAME.append(cl[0])
                self.Hindx[self.HNAME[-1]]=len(self.HNAME)-1
                self.Hindx[self.HNAME[-1].upper()]=len(self.HNAME)-1

                
            if corr==0 and cl[1]=='V' and cl[8]=="OK":
                self.VCO.append(float(cl[3])*1e-3)
                self.VCOerr.append(0.0)
                self.VNAME.append(cl[0])
                self.Vindx[self.VNAME[-1]]=len(self.VNAME)-1
                self.Vindx[self.VNAME[-1].upper()]=len(self.VNAME)-1

            if corr==1 and cl[1]=='H':
                self.Hcorr.append(cl[4])
                self.HcorrNAME.append(cl[0])
                self.Hindx[self.HcorrNAME[-1]]=len(self.HcorrNAME)-1
                self.Hindx[self.HcorrNAME[-1].upper()]=len(self.HcorrNAME)-1
            
            if corr==1 and cl[1]=='V':
                self.Vcorr.append(cl[4])
                self.VcorrNAME.append(cl[0])
                self.Vindx[self.VcorrNAME[-1]]=len(self.VcorrNAME)-1
                self.Vindx[self.VcorrNAME[-1].upper()]=len(self.VcorrNAME)-1 
        self.HCO=array(self.HCO)
        self.VCO=array(self.VCO)


    def rms(self):
      self.rmsx=std(self.HCO)
      self.rmsy=std(self.VCO)

    def minus(self, b):
      c=YASPtwiss()
      for bpm in self.VNAME:
        aindx=self.Vindx[bpm]
        
        try:
          bindx=b.Vindx[bpm]
          c.VNAME.append(bpm)
          c.Vindx[bpm]=len(c.VNAME)-1
          c.VCO.append(  self.VCO[aindx] -  b.VCO[bindx] )
          c.VCOerr.append( sqrt(self.VCOerr[aindx]**2 + b.VCOerr[bindx]**2) )
        except:
          print "Doing YASP V minus:", bpm, "is not in all files"
      
      for bpm in self.HNAME:
        aindx=self.Hindx[bpm]
        try:
          bindx=b.Hindx[bpm]
          c.HNAME.append(bpm)
          c.Hindx[bpm]=len(c.HNAME)-1
          c.HCO.append(  self.HCO[aindx] -  b.HCO[bindx] )
          c.HCOerr.append( sqrt( self.HCOerr[aindx]**2 + b.HCOerr[bindx]**2) ) 
          
        except:                                                               
          print "Doing YASP H minus:", bpm, "is not in all files"

      return c

    def AddModel(self, m1, m2):
      d=YASPtwiss()
      d.HS=[]
      d.VS=[]
      for i in range(len(self.HNAME)):
        bpm=self.HNAME[i]
        if "B1" in bpm:
            m=m1
        else:
            m=m2
        
        try:
          mindx=m.indx[bpm]
          flag=1
        except:
          print "Adding Model H:", bpm, "Not found in model"
          flag=0
          
        if flag==1:
          co=self.HCO[i]
          coerr=self.HCOerr[i]
          d.HNAME.append(bpm)
          d.Hindx[bpm]=len(d.HNAME)-1
          d.HCO.append(co)
          d.HCOerr.append(coerr)
          d.HS.append(m.S[mindx])
          
      for i in range(len(self.VNAME)):
        bpm=self.VNAME[i]
        if "B1" in bpm:
            m=m1
        else:
            m=m2
        try:
          mindx=m.indx[bpm]
          flag=1
        except:
          print "Adding Model V:", bpm, "Not found in model"
          flag=0 

        if flag==1:
          co=self.VCO[i]
          coerr=self.VCOerr[i]
          d.VNAME.append(bpm)
          d.Vindx[bpm]=len(d.VNAME)-1
          d.VCO.append(co)
          d.VCOerr.append(coerr)
          d.VS.append(m.S[mindx])

      return d



    def YASPwrite(self,ofile):
      try:
        self.HS[0]
      except:
        print "S is not defined, setting it to zero"
        self.HS=zeros(len(self.HNAME))
        self.VS=zeros(len(self.VNAME))
          
      f=open(ofile+'x',"w")
      print >>f, "* NAME    S    POS   ERR"
      print >>f, "$  %s    %le    %le   %le"
      for i in range(len(self.HNAME)):
        print >>f , self.HNAME[i], self.HS[i], self.HCO[i], self.HCOerr[i]
      f.close()
      f=open(ofile+'y',"w")
      print >>f, "* NAME    S    POS   ERR"
      print >>f, "$  %s    %le    %le   %le"
      for i in range(len(self.VNAME)):
        print >>f , self.VNAME[i], self.VS[i], self.VCO[i], self.VCOerr[i]
      f.close()

#########################################################
      # End of YASPtwiss
############################################################


def ModelDiff(a,b):
    s=[]
    x=[]
    y=[]
    for i in range(len(a.NAME)):
        s.append(a.S[i])
        x.append(-(a.X[i]-b.X[i])*1000)
        y.append(-(a.Y[i]-b.Y[i])*1000)

    return s,x,y


def ModelWrite(s,x,y,file,si=0,sf=999999):
  f=open(file,'w')
  for i in range(len(s)):
    if s[i] > si and s[i] < sf:
      print >>f, s[i],x[i],y[i]
  f.close()
  return 1
  




def  YASPaverage(list):
    "Average YASP classes"
    if len(list)==0:
        print "Nothing to average"
        return
    elif len(list)==1:
        return list[0]
    else:
        a=YASPtwiss()
        for bpm in list[0].VNAME:
            
            average=0
            rms=0
            count=0
            
            for files in list:
                flag=1
                
                try:
                    co=files.VCO[files.Vindx[bpm]]            
                    average=average+co
                    rms=rms+co**2
                    count=count+1
                except:
                    print "Doing V average:", bpm, "not found in all files"
                    flag=0

            if flag==1 and count > 0:
                
                a.VNAME.append(bpm)
                a.Vindx[bpm]=len(a.VNAME)-1
                a.VCO.append(average/count)
                a.VCOerr.append( sqrt(rms/count - average**2/count**2 + 1e-18) )

                
             
            count=0
            rms=0
            average=0
            flag=1
            for files in list:
                try:
                    co=files.HCO[files.Hindx[bpm]]            
                    average=average + co
                    rms=rms + co**2
                    count=count+1
                except:
                    print "Doing H average:", bpm, "not found in all files"
                    flag=0


            if flag==1 and count > 0:
                
                a.HNAME.append(bpm)
                a.Hindx[bpm]=len(a.HNAME)-1
                a.HCO.append(average/count)
                a.HCOerr.append( sqrt(rms/count - average**2/count**2  + 1e-18) )

    return a           




def GnuplotWrite(file,name,si,sf):
  
  f=open(file+".gplt",'w')
  
  print >>f, "set terminal postscript enhanced color solid 24"
  print >>f, "set output \""+name+".eps\""
  print >>f, "set multiplot"
  print >>f, "set size 1,0.54"
  print >>f, "set origin 0,0"
  print >>f, "set xrange[",si,":",sf,"]"
  print >>f, "set yrange [-3:3]"
  print >>f, "unset key"

  print >>f, "set rmargin 3"
  print >>f, "set lmargin 6"

  print >>f, "set xlabel \"Longitudinal location [m]\"" 
  print >>f, "set ylabel \"y [mm]\""

  print >>f, "p  \""+name+"y\" u 2:3:4 w e lt 3 lw 2, \"\" u 2:3 w l lt 3 lw 2, \""+name+".model\" u 1:3 w l lt 1 lw 2"

  print >>f, "set origin 0,0.48"
  print >>f, "set size 1,0.49"
  print >>f, "set xlabel \"\""
  print >>f, "set ylabel \"x [mm]\""
  print >>f, "set label \""+name+"\" at graph 0.8,1.06"
  print >>f, "p  \""+name+"x\" u 2:3:4 w e lt 3 lw 2, \"\" u 2:3 w l lt 3 lw 2, \""+name+".model\" u 1:2 w l lt 1 lw 2"
  f.close()
  


  













