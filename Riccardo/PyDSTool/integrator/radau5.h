#ifndef __RADAU5__
#define __RADAU5__

/* Calling sequence: 
   phaseDim, vfieldfunc, 
   tinit, y0, tend, h,   
   rtol, atol, itol, 
   vfieldjac, 
   ijac, mljac, mujac,
   vfieldmas,
   imas, mlmas, mumas,
   radau_solout, 
   iout,
   workArray, lenWorkArray, iWorkArray, leniWorkArray, 
   rpar, ipar, idid */

void radau5_(int*, void(*)(int*,double*,double*,double*,double*,int*),
	     double*, double*, double*, double*,
	     double*,double*,int*,
	     void(*)(int*,double*,double*,double*,int*,double*,int*),
	     int*,int*,int*,
	     void(*)(int*,double*,int*,double*,int*,double*,double*),
	     int*,int*,int*,
	     void(*)(int*,double*,double*,double*,double*,int*,int*,double*,
		     int*,int*),
	     int*,
	     double*,int*,int*,int*,
	     double*,int*,int*);

/* Calling sequence: 
   index + 1, time, contdata, contdim */

double contr5_(int*,double*,double*,int*);

/* This function is unused */
void dontr5_(double*,double*,int*,double*,double*,int*);

#endif
