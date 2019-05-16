import java.lang.Math;
import java.io.*;
import java.util.Locale;
//import org.knowm.xchart.XYChart;
//import org.knowm.xchart.QuickChart;
//import org.knowm.xchart.SwingWrapper;

public class Frisbeev2 
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
		double c=45; 	                            //maximum angle of attack
		double b=0.38*(c-a);
		double beta = b;
		
	double distancefinal = calculateDistance(beta);
	System.out.println("the final distance is: " + distancefinal);
		
	
        //double betaoptimal = goldenSearch(a,b,c); 
	    //System.out.println("the optimal angle of attack is " + betaoptimal);
	
	
	
	}


	public static double calculateLift(double angle, double beta)
	{
		double clift;
		double betarad =Math.toRadians(beta);
		double anglerad = Math.toRadians(angle);
		double angleattack = betarad - anglerad;
		
		clift = CL0 + CLalpha*angleattack;
		
		return clift;
	
	}
	public static double calculateLiftForce(double clift, double v1)
	{
		double forcelift;
		
		forcelift =0.5*rho*Math.pow(v1,2)*area*clift;
		
		return forcelift;
	}
	public static double calculateDrag(double angle, double beta)
	{
		double cdrag;
		double betarad = Math.toRadians(beta);
		double anglerad = Math.toRadians(angle);
		double alpha0rad= Math.toRadians(alpha0);
		double angleattack = betarad - anglerad;
		
		cdrag = CD0 + CDalpha*Math.pow((angleattack-alpha0rad),2);
		
		return cdrag;
	}
	public static double calculateDragForce(double cdrag, double v1)
	{
		double forcedrag;
		
		forcedrag = 0.5*cdrag*rho*area*Math.pow(v1,2);
		
		return forcedrag;
	}
		
	public static double goldenSearch(double a,double b,double c)
	{
		while (Math.abs(a-c)>(2.*a* 1.11e-16))
		{
			 double ya= calculateDistance(a);
			 double yb= calculateDistance(b);
			 double yc=calculateDistance(c);
			 double d; 
			
			
			if ((c-b) < (b-a))
			{
				 d = b- 0.38197*(b-a);
				double yd= -calculateDistance(d);
				
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
			double yd= -calculateDistance(d);
        
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
        System.out.println("the angle in degrees is"+ b);
		return b;
        
		
	}
	public static double accelerationX(double angle, double v1, double beta)
	{
        double angleradDrag = Math.toRadians(angle);
        double angleradLift = Math.toRadians(angle);
        double cdrag=calculateDrag(angle, beta);
        System.out.println("the coefficient of drag is: " + cdrag);
        double clift=calculateLift(angle,beta);
        System.out.println("the coefficient of lift is: " + clift);
		
		double aix = -Math.cos(angleradDrag)*calculateDragForce(cdrag, v1) - Math.cos(angleradLift)*calculateLiftForce(clift, v1);
		System.out.println("the acceleration in x is: " + aix);
		System.out.println("the drag force is: " + calculateDragForce(cdrag,v1));
		System.out.println("the lift force is: " + calculateLiftForce(clift,v1));
		
		
		return aix;
	}
	public static double accelerationY(double angle, double v1, double beta)
    {
        double angleradDrag = Math.toRadians(angle);
        double angleradLift = Math.toRadians(angle);
        double cdrag=calculateDrag(angle,beta);
        double clift=calculateLift(angle,beta);
		
		double aiy =  Math.sin(angleradDrag)*calculateDragForce(cdrag,v1) + Math.sin(angleradLift)*calculateLiftForce(clift,v1)-(m*g);
		System.out.println("the acceleration in y is: " + aiy);
		
		return aiy;
	}
		
	public static double calculateDistance(double beta)
	{
		double anglerad= Math.toRadians(beta);
				
		int imax= (int)(Tmax/dt);					//Maximum index
	
		double x[]= new double[imax];				//x-position of the frisbee
		double y[]= new double[imax];				//y-position of the frisbee
		
		x[0]=0.;								    //initial position of x in meter 
		y[0]=1.;								    //initial position of y in meter
		double yf = 0.;									//final position of y in meter
		System.out.println(x[0] + "     " + y[0]);

		double[] vy=new double[imax];				//velocity in the y-direction
		double[] vx=new double[imax];
		
		vx[0] = vi*Math.cos(anglerad);
		vy[0] = vi*Math.sin(anglerad);
		System.out.println("the initial velocity in x" + vx[0] + "the initial velocity in y" + vy[0]);
		
		double v1 = Math.sqrt(Math.pow(vx[0],2)+Math.pow(vy[0],2));
		
		double[] ax = new double[imax];
		double[] ay = new double[imax];
		
		
		ax[0]= accelerationX(beta,v1,beta);
		ay[0]= accelerationY(beta,v1,beta);
		System.out.println("the initial acceleration in x" + ax[0] + "the initial acceleration in y" + ay[0]);
		
		double[] angleupdate = new double[imax];
		angleupdate[0] = beta ;
		
		
		double deltay = yf - y[0];					//final y-position minus initial y position
		
		
		double distance = 0;
		
		//calculate the first value with Euler's method
			y[1] = y[0] + vy[0]*dt;
			x[1] = x[0] + vx[0];
	    	System.out.println("the postion at time 1 in x " + x[1] + "the position at time 1 in y " + y[1]);
			
			vx[1] = vx[0] + ax[0]*dt;
			vy[1] = vy[0] + ay[0]*dt;
			System.out.println("the velocity at time 1 in x " + vx[1] + "the velocity at time 1 in y " + vy[1]);
			
			v1 = Math.sqrt(Math.pow(vx[0],2)+Math.pow(vy[0],2));
			
			double angleinrad = Math.atan(vy[1]/vx[1]);
			angleupdate[1] = Math.toDegrees(angleinrad);
			System.out.println("The angle of the velocity is: " + angleupdate[1]);
			
			ax[1]= accelerationX(angleupdate[1],v1,beta);
		    ay[1]= accelerationY(angleupdate[1],v1,beta);
		    System.out.println("the initial acceleration in x " + ax[1] + "the initial acceleration in y " + ay[1]);
			
		int i=2;	
		
		while(y[i-1]>0)
		{
			double angleradian =  Math.atan(vx[i-1]/vy[i-1]);
			angleupdate[i] = Math.toDegrees(angleradian);
					
			ax[i] = accelerationX(angleupdate[i],v1,beta);
			ay[i] = accelerationY(angleupdate[i],v1,beta);
			System.out.println("The updated acceleration in x is: " + ax[i] + "The updated acceleration in y is: " + ay[i]);
			
			vx[i] = vx[i-1] + ax[i-1]*dt;
			vy[i] = vy[i-1] + ay[i-1]*dt;
			System.out.println("The updated velocity in x is: " + vx[i] + "the updated velocity in y is: " + vy[i]);
			
			v1 = Math.sqrt(Math.pow(vx[i],2)+Math.pow(vy[i],2));
			System.out.println("the v1 is: " + v1);
			
			x[i] = (2*x[i-1]) - x[i-2] + accelerationX(angleupdate[i],v1, beta)*Math.pow(dt,2);
			y[i] = (2*y[i-1]) - y[i-2] + accelerationY(angleupdate[i],v1, beta)*Math.pow(dt,2);
			
			x[i] = x[i-1] + vx[i-1]*dt;
			y[i] = y[i-1] + vy[i-1]*dt;
			System.out.println("the position in y is: " + y[i] + "the position in x is: " + x[i]);
			
			distance=x[i];
			i++;
		
		}
		System.out.println("distance is: " + distance);
		return distance;
        //System.out.println("value of the x position"+ x[]);
        //System.out.println("value of the y position"+ y[]);
	
	}
}