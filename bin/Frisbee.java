import java.lang.Math;
import java.io.*;
import java.util.Locale;
//import org.knowm.xchart.XYChart;
//import org.knowm.xchart.QuickChart;
//import org.knowm.xchart.SwingWrapper;

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
	
	public static final double vi=14.;			//initial velocity of average throw(m/s)
	public static void main(String[] args) 
	{	
	
	
		double a=0;                                 //minimum angle of attack
		double c=90; 	                            //maximum angle of attack
		double b=0.38*(c-a);
		
		
	
		double answer = goldenSearch(a,b,c); 
		System.out.println(answer);
		
	/////////////GRAPH/////////////////
	
		
	//	XYChart chart = QuickChart.getChart("Distance in x relative to the height of the frisbee","X","Y","y(x)",x,y);
	//	new SwingWrapper(chart).displayChart();
		
		
	
	
	}


	public static double calculateLift(double angle)
	{
		double clift;
		
		clift = CL0 + CLalpha*angle;
		
		return clift;
	
	}
	public static double calculateLiftForce(double clift, double vi)
	{
		double forcelift;
		
		forcelift = 0.5*rho*Math.pow(vi,2)*area*clift;
		
		return forcelift;
	}
	public static double calculateDrag(double angle)
	{
		double cdrag;
		
		cdrag = CD0 + CDalpha*Math.pow((angle-alpha0),2);
		
		return cdrag;
	}
	public static double calculateDragForce(double cdrag, double vi)
	{
		double forcedrag;
		
		forcedrag = -0.5*cdrag*rho*area*Math.pow(vi,2);
		
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
	public static double accelerationX(double angle)
	{
        double angleradDrag = Math.toRadians(angle);
        double angleradLift = Math.toRadians(180-90-angle);
        double cdrag=calculateDrag(angle);
        double clift=calculateLift(angle);
		
		double aix = Math.cos(angleradDrag)*calculateDragForce(cdrag, vi) + Math.cos(angleradLift)*calculateLiftForce(clift, vi);
		
		
		return aix;
	}
	public static double accelerationY(double angle)
    {
        double angleradDrag = Math.toRadians(angle);
        double angleradLift = Math.toRadians(180-90-angle);
        double cdrag=calculateDrag(angle);
        double clift=calculateLift(angle);
		
		double aiy = Math.sin(angleradDrag)*calculateDragForce(cdrag,vi) + Math.sin(angleradLift)*calculateLiftForce(clift,vi)-(m*g);
		
		return aiy;
	}
		
	public static double calculateDistance(double angle)
	{
		double anglerad= Math.toRadians(angle);
				
		int imax= (int)(Tmax/dt);					//Maximum index
	
		double x[]= new double[imax];				//x-position of the frisbee
		double y[]= new double[imax];				//y-position of the frisbee
		
		x[0]=0.;								    //initial position of x in meter 
		y[0]=1.;								    //initial position of y in meter
		double yf = 0.;									//final position of y in meter

		double[] vy=new double[imax];				//velocity in the y-direction
		double[] vx=new double[imax];
		
		vx[0] = vi*Math.cos(anglerad);
		vy[0] = vi*Math.sin(anglerad);
		
		double[] ax = new double[imax];
		double[] ay = new double[imax];
		
		ax[0]= accelerationX(angle);
		ay[0]= accelerationY(angle);
		
		double[] angleupdate = new double[imax];
		angleupdate[0] = angle ;
		
		
		double deltay = yf - y[0];					//final y-position minus initial y position
		
		
		double distance = 0;
		
		//calculate the first value with Euler's method
			y[1] = y[0] + vy[0]*dt;
			x[1] = x[0] + vx[0];
			
			vx[1] = vx[0] + accelerationX(angle)*dt;
			vy[1] = vy[0] + accelerationY(angle)*dt;
			
			double angleinrad = Math.atan(vy[1]/vx[1]);
			angleupdate[1] = Math.toDegrees(angleinrad);
		
		for(int i=2; y[i]> yf ; i++)
		{
			double angleradian =  Math.atan(vx[i-1]/vy[i-1]);
			angleupdate[i] = Math.toDegrees(angleradian);
					
			ax[i] = accelerationX(angleupdate[i]);
			ay[i] = accelerationY(angleupdate[i]);
			
			vx[i] = (1 / (2*m))  * rho * Math.pow(vx[i-1],2) * area * calculateDrag(angleupdate[i]) * dt;
			vy[i] = (g + (1/(2*m)))*rho*Math.pow(vx[i-1],2)*area*calculateLift(angleupdate[i])*dt;
			
			x[i] = 2*x[i-1] - x[i-2] + accelerationX(angleupdate[i])*Math.pow(dt,2);
			y[i] = 2*y[i-1] - y[i-2] + accelerationY(angleupdate[i])*Math.pow(dt,2);
			
			distance=x[i];
		
		}
		System.out.println("distance is: " + distance);
		return distance;
	
	}

}	