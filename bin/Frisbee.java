import java.lang.Math;
import java.io.*;
import java.util.Locale;
import ptolemy.plot.*;

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
	public static final double g=9.81;			//acceleration due to gravity (m/s^2)
	public static final double m=0.175; 		//mass of a frisbee (in kg)
	
	
	
	double a=0;                                 //minimum angle of attack
	double c=90; 	                            //maximum angle of attack
	double b=0.38*(c-a);
	
	
	int imax= (int)(Tmax/dt);					//Maximum index
	
	double[] x= new double[imax];				//x-position of the frisbee
	double[] y= new double[imax];				//y-position of the frisbee

	double[] vy=new double[imax];				//velocity in the y-direction
	double[] vx=new double[imax];
	
	x[0]=0.;								//initial position of x in meter 

	y[0]=1.;								//initial position of y in meter
	y[imax]=0.;								//final position of y in meter
	
	
		
	double vi=14.;								//initial velocity of average throw(m/s)
	public static void main(String[] args)
	{	
		double answer = goldenSearch(a,b,c); 
		System.out.println(answer);
		
		Plot plotObj = new Plot (); // Create Plot object
		plotObj.setTitle ("Distance in y vs Distance in x");
		plotObj.setXLabel ("Distance in x");
		plotObj.setYLabel ("Distance in y");
	// plotObj.setSize (400, 300) ;
	// plotObj.setXRange(Xmin, Xmax );
	// plotObj.addPoint(int Set, double x, double y, boolean connect)
		double xS tep = (Xmax−Xmin) / Npoint ;
	// Plotting loop
		for (double x=Xmin ; x<=Xmax ; x+=xS tep) 
		{
			double y = Math.sin(x)*Math.sin(x);
			plotObj.addPoint (0,x,y,true);
		}
		
    PlotApplication app = new PlotApplication (plotObj) ; // Display
}


	}
<<<<<<< HEAD
	public static double calculateCoefficientLift(double angle,double vx)
=======
	public static double calculateLift(double angle,double v)
>>>>>>> c2c7a628ed83525cab71d89a6ffaba6db7c5709e
	{
		double clift;
		
		clift = CL0 + CLalpha*angle;
<<<<<<< HEAD
		
		return clift;
	
	}
	public static double calculateLiftForce(double clift, double vx)
	{
		double forcelift;
		
		forcelift = 0.5*rho*Math.pow(vx,2)*area*clift;
=======
		forcelift = 0.5*rho*Math.pow(v,2)*area*clift;
>>>>>>> c2c7a628ed83525cab71d89a6ffaba6db7c5709e
		
		return forcelift;
	}
<<<<<<< HEAD
	public static double calculateCoefficientDrag(double angle,double vx)
=======
	public static double calculateDrag(angle,v)
>>>>>>> c2c7a628ed83525cab71d89a6ffaba6db7c5709e
	{
		double cdrag;
		
		cdrag = CD0 + CDalpha*Math.pow((angle-alpha0),2);
<<<<<<< HEAD
		
		return cdrag;
	}
	public static double calculateDragForce(double cdrag, double vx)
	{
		double forcedrag;
		
		forcedrag = -0.5*cdrag*rho*area*Math.pow(vx,2);
=======
		forcedrag = -0.5*cdrag*rho*area*Math.pow(v,2);
>>>>>>> c2c7a628ed83525cab71d89a6ffaba6db7c5709e
		
		return forcedrag;
	}
		
	public static double goldenSearch(double a,double b,double c)
	{
		while (Math.abs(a-c)>(2.*a* 1.11e-16))
		{
			 double ya= calculateDistance(a);
			 double yb= calculateDistance(b);
			 double yc= calculateDistance(c);
			 double d; 
			
			
			if ((c-b) < (b-a))
			{
				 d = b- 0.38197*(b-a);
				double yd= calculateDistance(d);
				
				if (yd<yb)
				{
				c=b;
				b=d;
				}
				else
				{
                a=d;
				}
				System.out.println(b);
            }
			else 
			{
			d = b+ 0.38197*(c-b);
			double yd= calculateDistance(d);
        
				if (yd<yb)
				{
				a = b;
				b = d;
				}
				else
				{
                c = d;
				}
				System.out.println(b);
			}
			
		}
		return b;
		
	}
	public static double calculateDistance(double angle)
	{
		double anglerad= Math.toRadians(angle);
		vx[0] = vi*Math.cos(anglerad);
		vy[0] = vi*Math.sin(anglerad);
		double deltay= y[imax] - y[0];
		double[] i = new  double[imax];
		double cdrag = calculateCoefficientDrag(angle, vx);
		double clift = calculateCoefficientLift(angle,vx);
		
		double aix = calculateDrag(angle) / m;
		double aiy = (calculateLift(angle)-calculateDrag(angle))/m;
		
		for(int i=1; y[i]> y[imax] ; i++)
		{
			//calculate the first value with Euler's method
			y[i] = y[i-1] + vy[i-1]*dt;
			
			vx[i] = (1 / (2*m))  * rho * Math.pow(vx[i-1],2) * area * cdrag * dt;
			vy[i] = (g + (1/(2*m))*rho*Math.pow(vx[i-1],2)*area*clift)*dt;
			
			
			t[i] = (-vy[i] + Math.sqrt(Math.pow(vy[i],2)-4*(0.5*-g)*deltay))/ (2*deltay);
		
			if(t<0)
			{	
			t[i] = (-vy[i] - Math.sqrt(Math.pow(vy[i],2)-4*(0.5*-g)*deltay))/ (2*deltay);
			}
		
			//calculate the first value with Euler's method
		
			x[i] = x[i-1] + vx[i-1];
			double distance=x[i];
		}
		
		return distance;
	
	}

}	