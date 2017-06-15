/*
	Single Sided Gaseous Lid Driven Micro-flow using Lattice Boltzmann Method
	Author : Isac Rajan
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdio>
//#include <string>
//#include <sstream>

#define NX 301 // No of grids in X direction
#define NY 301 // No of grids in Y direction

using std::vector;
using std::cout;
using std::endl;
using std::ofstream;
//using std::string;
//using std::ostringstream;

typedef double mydoub;
typedef int myint;
typedef vector<vector<vector<double> > > vdouble3;
typedef vector<vector<double> > vdouble2;

myint TSTEPS = 100000; //Max number of timesteps
myint cx[9] = {0,1,0,-1,0,1,-1,-1,1};
myint cy[9] = {0,0,1,0,-1,1,1,-1,-1};
mydoub w[9] = {4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36}; //weights
mydoub rhoo = 1.0;
mydoub Kn = 0.01; // Knudsen Number
mydoub tau = Kn*(NY-1); // Relaxation Time
mydoub omega = 1.0/tau;
mydoub sigma = 0.7; // TMAC Tangetial Momentum Accomodation Coefficient
mydoub uo = 0.001;
vdouble3 f, feq;
vdouble2 u, v, rho, u_old, v_old;

void Initialize();
void AssignOldVector();
void Collision();
void Stream();
void BoundaryCondition();
void ComputeField();
mydoub ComputeError();
void InitializeLid();
void WriteToFile();

int main()
{
	int t;
	clock_t ct; // Variable to measure CPU time of execution
	ct = clock();
	mydoub error;
	Initialize();

	//Main Loop Starts
	for(t=1; t<=TSTEPS; t++)
	{
		AssignOldVector();
		Collision();
		Stream();
		BoundaryCondition();
		ComputeField();
		error = ComputeError();
		InitializeLid();

		if(t%10==0)
		{
			cout << "Timestep = " << t << endl;
			cout << "Time taken till now = " << (double)((-ct + clock())/CLOCKS_PER_SEC)/60.0 << " mins" <<endl;
			cout << "Error = "<< error <<endl;
			if (error <= 1e-6)
				break;
		}

		if(t%500==0)
			WriteToFile();
	}
	cout << "Final Error = " << error << " at Timestep = " << t << "\n\n";
	WriteToFile();

	return 0;
}

void Initialize()
{
	//Initializing rho, u, v vectors
	rho.resize(NX);
	u.resize(NX); u_old.resize(NX);
	v.resize(NX); v_old.resize(NX);
	for(int i=0; i<NX; i++)
	{
		rho[i].resize(NY);
		u[i].resize(NY); u_old[i].resize(NY);
		v[i].resize(NY); v_old[i].resize(NY);
	}
	for(int i=0; i<NX; i++)
		for(int j=0; j<NY; j++)
			{
				rho[i][j] = rhoo;
				u[i][j] = 0.0; u_old[i][j] = 0.0;
				v[i][j] = 0.0; v_old[i][j] = 0.0;
			}

	//Initializing f and feq vector to zeroes
	f.resize(9);
	feq.resize(9);
	for(int k=0; k<9; k++)
	{
		f[k].resize(NX);
		feq[k].resize(NX);
	}
	for(int k=0; k<9; k++)
		for(int i=0; i<NX; i++)
		{
			f[k][i].resize(NY);
			feq[k][i].resize(NY);
		}
	for(int k=0; k<9; k++)
		for(int i=0; i<NX; i++)
			for(int j=0; j<NY; j++)
			{
				f[k][i][j] = 0.0;
				feq[k][i][j] = 0.0;
			}

	//Initializing the velocity of lid
	for (int i=0; i<NX; i++)
	{
		u[i][NY-1] = uo;
		v[i][NY-1] = 0.0;
	}

	//Initializing the f with feq at t=0
	for(int j=0; j<NY; j++)
	{
		for(int i=0; i<NX; i++)
		{
			mydoub t1;
			t1=u[i][j]*u[i][j]+v[i][j]*v[i][j];
			for(int k=0; k<=8; k++)
			{
				mydoub t2;
				t2=(u[i][j]*cx[k])+(v[i][j]*cy[k]);
				feq[k][i][j]=rho[i][j]*w[k]*(1.0+(3.0*t2)+(4.5*t2*t2)-(1.50*t1));
				f[k][i][j] = feq[k][i][j];
			}
		}
	}
}

void AssignOldVector()
{
	//Assignment of old arrays
	for(int j=0; j<NY; j++)
		for(int i=0; i<NX; i++)
		{
			u_old[i][j] = u[i][j];
			v_old[i][j] = v[i][j];
		}
}

void Collision()
{
	for(int j=0; j<NY; j++)
		for(int i=0; i<NX; i++)
		{
			mydoub t1;
			t1=u[i][j]*u[i][j]+v[i][j]*v[i][j];
			for(int k=0; k<=8; k++)
			{
				mydoub t2;
				t2=u[i][j]*cx[k]+v[i][j]*cy[k];
				feq[k][i][j]=rho[i][j]*w[k]*(1.0+(3.0*t2)+(4.5*t2*t2)-(1.50*t1));
				f[k][i][j]=omega*feq[k][i][j]+(1.0-omega)*f[k][i][j];
			}
		}
}

void Stream()
{
    /*
    for(int j=0; j<NY; j++)
    {
        for(int i=NX-1; i>=1; i--) //Right to Left
            f[1][i][j]=f[1][i-1][j];

        for(int i=0; i<NX-1; i++) //Left to Right
            f[3][i][j]=f[3][i+1][j];
    }
    for(int j=NY-1; j>=1; j--) //Top to Bottom
    {
        for(int i=0; i<NX; i++)
            f[2][i][j]=f[2][i][j-1];

        for(int i=NX-1; i>=1; i--)
            f[5][i][j]=f[5][i-1][j-1];

        for(int i=0; i<NX-1; i++)
            f[6][i][j]=f[6][i+1][j-1];
    }
    for(int j=0; j<NY-1; j++) //Bottom to Top
    {
        for(int i=0; i<NX; i++)
            f[4][i][j]=f[4][i][j+1];

        for(int i=0; i<NX-1; i++)
            f[7][i][j]=f[7][i+1][j+1];

        for(int i=NX-1; i>=1; i--)
            f[8][i][j]=f[8][i-1][j+1];
    }
    */
    int i,j;
     for(i=0;i<NX;i++)
      for(j=0;j<(NY-1);j++)
        f[0][i][j]=f[0][i][j];

     for(i=NX-1;i>0;i--)
      for(j=0;j<(NY-1);j++)
       f[1][i][j]=f[1][i-1][j];

     for(i=0;i<NX;i++)
      for(j=NY-2;j>0;j--)
       f[2][i][j]=f[2][i][j-1];

     for(i=0;i<(NX-1);i++)
      for(j=0;j<(NY-1);j++)
       f[3][i][j]=f[3][i+1][j];

     for(i=0;i<NX;i++)
      for(j=0;j<(NY-1);j++)
       f[4][i][j]=f[4][i][j+1];

     for(i=NX-1;i>0;i--)
      for(j=NY-2;j>0;j--)
       f[5][i][j]=f[5][i-1][j-1];

     for(i=0;i<(NX-1);i++)
      for(j=NY-2;j>0;j--)
       f[6][i][j]=f[6][i+1][j-1];

     for(i=0;i<(NX-1);i++)
      for(j=0;j<(NY-1);j++)
       f[7][i][j]=f[7][i+1][j+1];

     for(i=NX-1;i>0;i--)
      for(j=0;j<(NY-1);j++)
       f[8][i][j]=f[8][i-1][j+1];

}

void BoundaryCondition()
{
	mydoub rhon;
	for(int j=1; j<NY-1; j++)
	{
		//Bounce Back on West Boundary
		f[1][0][j]=f[3][0][j];
		f[5][0][j]=sigma*f[7][0][j] + (1.0-sigma)*f[6][0][j];
		f[8][0][j]=sigma*f[6][0][j] + (1.0-sigma)*f[7][0][j];

		//Bounce Back on East Boundary
		f[3][NX-1][j]=f[1][NX-1][j];
		f[7][NX-1][j]=sigma*f[5][NX-1][j] + (1.0-sigma)*f[8][NX-1][j];
		f[6][NX-1][j]=sigma*f[8][NX-1][j] + (1.0-sigma)*f[5][NX-1][j];
	}
	//Bounce Back on South Boundary
	for(int i=1; i<NX-1; i++)
	{

		f[2][i][0]=f[4][i][0];
		f[5][i][0]=sigma*f[7][i][0] + (1.0-sigma)*f[8][i][0];
		f[6][i][0]=sigma*f[8][i][0] + (1.0-sigma)*f[7][i][0];
	}
	//Moving Lid, North boundary
	for(int i=0; i<NX; i++)
	{
		for(int k=0; k<=8; k++)
		{
			f[k][i][NY-1] = feq[k][i][NY-1]; //setting the lid distribution to eqb
		}
	}
	for(int k=0; k<=8; k++) //one corner of the lid
	{
		f[k][0][0] = feq[k][0][0];
	}

	for(int k=0; k<=8; k++) //another corner of the lid
	{
		f[k][NX-1][0] = feq[k][NX-1][0];
	}
}

void ComputeField()
{
	mydoub ssum, usum, vsum;
	//Computation of the density
	for(int j=0; j<NY; j++)
	{
		for(int i=0; i<NX; i++)
		{
			ssum=0.0;
			usum=0.0;
			vsum=0.0;
			for(int k=0; k<=8; k++)
			{
				ssum+=f[k][i][j];
				usum+=f[k][i][j]*cx[k];
				vsum+=f[k][i][j]*cy[k];
			}
			rho[i][j]=ssum;
			u[i][j]=usum/rho[i][j];
			v[i][j]=vsum/rho[i][j];
		}
	}
}

mydoub ComputeError()
{
	mydoub err1, err2, err;
	err=err1=err2=0.0;
	mydoub ur[NX][NY];

	for(int i=0; i<NX; i++)
	{
		for(int j=0; j<NY; j++)
		{
			ur[i][j] = u[i][j]*u[i][j] + v[i][j]*v[i][j];
			err1 += (u[i][j]-u_old[i][j])*(u[i][j]-u_old[i][j]) + (v[i][j]-v_old[i][j])*(v[i][j]-v_old[i][j]);
			err2 += ur[i][j];
		}
	}
	err = sqrt(err1)/sqrt(err2);

	return err;
}

void InitializeLid()
{
	//Initializing the velocity of the top lid
	for (int i=0; i<NX; i++)
	{
		u[i][NY-1]=uo;
		v[i][NY-1]=0.0;
	}
}

void WriteToFile()
{
    /*
    ostringstream k, sig;
    k << Kn;
    sig << sigma;
	string str1 = "XYUV_Kn" + k.str() + "sig" + sig.str() + ".dat";
	string str2 = "YcUc_Kn" + k.str() + "sig" + sig.str() + ".dat";
	string str3 = "XcVc_Kn" + k.str() + "sig" + sig.str() + ".dat";
	*/

	//Writing data to file, X,Y,U,V -----------------------------
	//To be postprocessed using Tecplot
	ofstream myfile ("XYUV.dat");
	myfile << "VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n" << endl;
	myfile << "ZONE  F=POINT\n" << endl;
	myfile << "I=" << NX << ", J=" << NY << endl;
	double xpos, ypos, Dx, Dy;
	Dx=1.0/NX;
	Dy=1.0/NY;
	for(int j=0; j<NY; j++)
	{
		for(int i=0; i<NX; i++)
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
	for(int j=0; j<NY; j++)
	{
		myfile2 << j*Dy << "\t" << (u[(NX-1)/2][j])/(2*uo) << endl;
	}
	myfile2.close();

	//Writing Xc vs Vc data to file
	ofstream myfile3 ("XcVc.dat");
	for(int i=0; i<NX; i++)
	{
		myfile3 << i*Dx << "\t" << (v[i][(NY-1)/2])/(2*uo) << endl;
	}
	myfile3.close();
}
