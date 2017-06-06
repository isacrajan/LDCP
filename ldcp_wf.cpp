#include <iostream>
#include <ctime>
#include <fstream>

#define N 290 //No of grids in X-direction
#define M 290 //No of grids in Y-direction

using namespace std;

int main()
{
	clock_t ct;
	ct=clock(); //Variable to measure CPU time

	int tsteps, t;
	tsteps=10000;	//no of iterations

	double f[9][N+1][M+1], feq[9][N+1][M+1];
	double rho[N+1][M+1], u[N+1][M+1], v[N+1][M+1];
	
	//Lattice Weights
	double w[9]={4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36};

	//lattice velocities for D2Q9
	int cx[9]={0,1,0,-1,0,1,-1,-1,1}, cy[9]={0,0,1,0,-1,1,1,-1,-1};

	double t1,t2;
	double uo=0.1, Re=100; //lattice characteristic velocity and lattice Reynolds number
	double rhoo=1.0;
	double dx=1.0, dy, dt=1.0;
	dy=dx;
	double nu; //Kinematic viscosity
	nu=(double)uo*M/Re;
	cout << "The Reynolds Number is " << Re << "." << endl;

	double omega; // nu = (2tau-1)/6 * dx*dx/3*dt and omega=1/tau
	omega=2.0/(1.0 + 6.0*(nu*dt/(dx*dx)));
	cout << "The value of tau is " << 1.0/omega << " and omega is " << omega << endl;

	//Initializing density and velocities
	for(int j=0; j<=M; j++)
	{
		for(int i=0; i<=N; i++)
		{
			rho[i][j]=rhoo;
			u[i][j]=0.0; 
			v[i][j]=0.0;
		}
	}

	//Initializing the velocity of the lid
	for (int i=0; i<=N; i++)
	{
		u[i][M]=uo;
		v[i][M]=0.0;
	}

	//Initializing the f with feq at t=0
	for(int j=0; j<=M; j++)
	{
			for(int i=0; i<=N; i++)
			{
				t1=u[i][j]*u[i][j]+v[i][j]*v[i][j];
				for(int k=0; k<=8; k++)
				{
					t2=u[i][j]*cx[k]+v[i][j]*cy[k];
					f[k][i][j]=rho[i][j]*w[k]*(1.0+3.0*t2+4.5*t2*t2-1.50*t1);
				}
			}
	}

	//Main Loop Starts
	for(t=1; t<=tsteps; t++)
	{		
		//COLLISION STEP
		for(int j=0; j<=M; j++)
			for(int i=0; i<=N; i++)
			{
				t1=u[i][j]*u[i][j]+v[i][j]*v[i][j];
				for(int k=0; k<=8; k++)
				{
					t2=u[i][j]*cx[k]+v[i][j]*cy[k];
					feq[k][i][j]=rho[i][j]*w[k]*(1.0+3.0*t2+4.5*t2*t2-1.50*t1);
					f[k][i][j]=omega*feq[k][i][j]+(1.0-omega)*f[k][i][j];
				}
			}

		//STREAMING STEP
		for(int j=0; j<=M; j++)
		{
			for(int i=N; i>=1; i--) //Right to Left
				f[1][i][j]=f[1][i-1][j];

			for(int i=0; i<=N-1; i++) //Left to Right
				f[3][i][j]=f[3][i+1][j];
		}
		for(int j=M; j>=1; j--) //Top to Bottom
		{
			for(int i=0; i<=N; i++)
				f[2][i][j]=f[2][i][j-1];

			for(int i=N; i>=1; i--)
				f[5][i][j]=f[5][i-1][j-1];

			for(int i=0; i<=N-1; i++)
				f[6][i][j]=f[6][i+1][j-1];
		}
		for(int j=0; j<=M-1; j++) //Bottom to Top
		{
			for(int i=0; i<=N; i++)
				f[4][i][j]=f[4][i][j+1];

			for(int i=0; i<=N-1; i++)
				f[7][i][j]=f[7][i+1][j+1];

			for(int i=N; i>=1; i--)
				f[8][i][j]=f[8][i-1][j+1];
		}

		//BOUNDARY CONDITIONS
		double rhon;
		for(int j=0; j<=M; j++)
		{
			//Bounce Back on West Boundary
			f[1][0][j]=f[3][0][j];
			f[5][0][j]=f[7][0][j];
			f[8][0][j]=f[6][0][j];

			//Bounce Back on East Boundary
			f[3][N][j]=f[1][N][j];
			f[7][N][j]=f[5][N][j];
			f[6][N][j]=f[8][N][j];
		}
		//Bounce Back on South Boundary
		for(int i=0; i<=N; i++)
		{

			f[2][i][0]=f[4][i][0];
			f[5][i][0]=f[7][i][0];
			f[6][i][0]=f[8][i][0];
		}
		//Moving Lid, North boundary
		for(int i=1; i<=N-1; i++)
		{
			rhon=f[0][i][M]+f[1][i][M]+f[3][i][M]+2.0*(f[2][i][M]+f[6][i][M]+f[5][i][M]);
			f[4][i][M]=f[2][i][M];
/* 			f[8][i][M]=f[6][i][M] + uo*rhon/6.0;
			f[7][i][M]=f[5][i][M] - uo*rhon/6.0; */
			f[7][i][M]=f[5][i][M] + 0.5*(f[1][i][M]-f[3][i][M]) - 0.5*rhon*uo;
			f[8][i][M]=f[6][i][M] + 0.5*(f[3][i][M]-f[1][i][M]) + 0.5*rhon*uo;
		}

		//COMPUTATION OF RHO, U, V
		double ssum, usum, vsum;
		//Computation of the density
		for(int j=0; j<=M; j++)
		{
			for(int i=0; i<=N; i++)
			{
				ssum=0.0;
				for(int k=0; k<=8; k++)
				{
					ssum+=f[k][i][j];
				}
				rho[i][j]=ssum;
			}
		}

/* 		for(int i=1; i<=N; i++)
			rho[i][M]=f[0][i][M]+f[1][i][M]+f[3][i][M]+2.0*(f[2][i][M]+f[6][i][M]+f[5][i][M]); */

		//Computation of u, v velocities
		for(int i=0; i<=N; i++)
		{
			for(int j=0; j<=M; j++)
			{
				usum=0.0;
				vsum=0.0;
				for(int k=0; k<=8; k++)
				{
					usum+=f[k][i][j]*cx[k];
					vsum+=f[k][i][j]*cy[k];
				}
				u[i][j]=usum/rho[i][j];
				v[i][j]=vsum/rho[i][j];
			}
		}
				
		//Initializing the velocity of the lid
		for (int i=1; i<=N-1; i++)
		{
			u[i][M]=uo;
			v[i][M]=0.0;
		}

		if(t%50==0)
			cout << "Timestep = " << t << endl;

	}
	//End of Main Loop

	//Writing data to file, X,Y,U,V
	//To be postprocessed using Tecplot
	ofstream myfile ("XYUV.dat");
	myfile << "VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n" << endl;
	myfile << "ZONE  F=POINT\n" << endl;
	myfile << "I=" << M+1 << ", J=" << N+1 << endl;
	double xpos, ypos, Dx, Dy;
	Dx=1.0/N;
	Dy=1.0/M;
	for(int j=0; j<=M; j++)
	{
		for(int i=0; i<=N; i++)
		{
			xpos=i*Dx;
			ypos=j*Dy;
			myfile << xpos << "\t" << ypos << "\t" << u[i][j]/uo << "\t" << v[i][j]/uo << endl;
		}
	}
	myfile.close();

	//Writing center line data to file, Yc vs Uc
	//Will be plotted using Python
	ofstream myfile2 ("YcUc.dat");
	for(int j=0; j<=M; j++)
	{
		myfile2 << j*Dy << "\t" << (u[M/2][j]+u[-1+M/2][j])/2 << endl;
	}
	myfile2.close();

	//Writing Xc vs Vc data to file
	ofstream myfile3 ("XcVc.dat");
	for(int i=0; i<=N; i++)
	{
		myfile3<< i*Dx << "\t" << (v[i][N/2] + v[i][-1+ N/2])/2 << endl;
	}
	myfile3.close();


	//Displaying CPU time
	ct=clock()-ct;
	cout << "\n\nTime Taken = " << (double)(ct/CLOCKS_PER_SEC)/60.00 << " mins." << endl;
	return 0;
}
