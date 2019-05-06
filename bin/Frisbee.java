importjava.lang.Math;
importjava.io.*;

public class Frisbee
{
	public static final double dt=0.001;		//time in seconds
	public static final double Tmax= 60.;		//maximum time in seconds
	public static final double rho=1.23;		//density of air at sea level in kg per meter cube
	public static final double area=0.0531;		//Area of the frisbee with a diameter of 0.26 meter in meter square
	public static final double CL0=0.15;		//y-intercept of the linear function when the angle of attack is 0
	public static final double CLalpha=1.4;		//slope of the graph of the linear function of angle of attack
	public static final double CD0=0.08;		//minimum drag due to friction and pressure drag
	public static final double CDalpha=2.72;	//induced drag(varies with lift)
	public static final double alpha0=-4.;		//angle of attack that produces the least lift(in degree)
	public static final double Fg=-1.71675;		//force due to gravity corresponds to a mass of 0.175 kg and an acceleration due to gravity of 9.81 m/s^2
	
	double vi=14.;								//initial velocity of average throw(m/s)
	public static void main(String[] args)
	{	
	int imax= (int)(Tmax/dt);					//Maximum index
	
	double[] x= new double[imax];				//x-position of the frisbee
	double[] y= new double[imax];				//y-position of the frisbee
	double[] vx= new double[imax];				//velocity in the x-direction
	double[] vy=new double[imax];				//velocity in the y-direction
	
	x[0]=0.;								//initial position of x in meter 
	y[0]=1.;								//initial position of y in meter
	vx[0]=0.;
	vy[0]=0.;
	
	}
	public static double calculateLift(angle,vx)
	{
		double clift;
		double forcelift;
		
		clift = CL0 + CLalpha*angle;
		forcelift = 0.5*rho*Math.pow(vx,2)*area*clift;
		
		return forcelift;
	}
	public static double calculateDrag()
	{
	}
	
	
	