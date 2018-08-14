#include "udf.h"
DEFINE_PROFILE(outlet_pressure, thread, position)
{
	real t;
	real x[ND_ND];
	real z0;
	real y0;
	face_t f;
	
	begin_f_loop(f, thread)
	{
		t=RP_Get_Real("flow-time");
		F_CENTROID(x,f,thread);		
		y0 = x[1];
		z0 = x[2];
		 if (y0*y0+z0*z0<=0.040*0.040)	
		  { 
			if (t>0 && t<=0.01)
				{
				F_PROFILE(f,thread,position)=100000+0.1*101325.0*exp(-4.0*200.0*200.0*log(10.0)*(t-0.005)*(t-0.005))*sin(8000.0*3.14159*(t-0.005));
				}
			else
				{
				F_PROFILE(f,thread,position)=100000;
				}
		   }
		else
		  {
				F_PROFILE(f,thread,position)=100000;			
		  } 
	}
	end_f_loop(f,thread)
}

/* DEFINE_PROPERTY(cell_density,cell,thread) 
	{

		real density;

		real pressure=C_P(cell,thread);


		density=998.2*(1+pressure/1e-9);

		return density;
	} */
	

#define BMODULUS 1.0e9
#define rho_ref 998.2
#define p_ref 101325
DEFINE_PROPERTY(compressiblefluid_density, c, t)
{ 
	real rho; 
	real p, dp;
	real p_operating;
	p_operating = RP_Get_Real ("operating-pressure");
	p= C_P(c,t)+p_operating;
	dp=p-p_ref;
	rho= rho_ref*(1.0+dp/BMODULUS);
	return rho;
}


DEFINE_PROPERTY(sound_speed, c,t)
{
    real a; 
    real p, dp,p_operating;
    p_operating = RP_Get_Real ("operating-pressure");

    p = C_P(c,t) + p_operating;
    dp = p-p_ref; 
    a = (1.-dp/BMODULUS)*sqrt(BMODULUS/rho_ref);   
    return a;
}