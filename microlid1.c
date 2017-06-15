#include<stdio.h>
#include<math.h>

#define IMAX 251

int its;
double U,tou,Kn=0.1,cs,visc,Re=0.00024;
double rho[IMAX][IMAX],u[IMAX][IMAX],v[IMAX][IMAX],un[IMAX][IMAX],vn[IMAX][IMAX],ur[IMAX][IMAX];
double dt,t,sig,f[9][IMAX][IMAX],feq[9][IMAX][IMAX],w[9],err,e[9][2];

initialize()
{
 int i,j,k;
 float temp;
 
sig=0.1;
 
 tou = Kn*(IMAX-1);
 cs=(1/sqrt(3.00));
 U= 0.001;

 for(i=0;i<IMAX;i++)
  for(j=0;j<IMAX;j++)
   {
    rho[i][j]=1.0;
    if(j!=(IMAX-1))
    u[i][j]=0;
    else
    {
     if(i==0|| i==IMAX-1)
     u[i][j]=0;
     else
     u[i][j]=U;
    }
     v[i][j]=0;
    ur[i][j]=(u[i][j]*u[i][j])+(v[i][j]*v[i][j]);
   }

 w[0]=0.444444;
 w[1]=w[2]=w[3]=w[4]=0.111111;
 w[5]=w[6]=w[7]=w[8]=0.027778;
 
 e[0][0]=e[0][1]=0;
 e[1][0]=1;e[1][1]=0;
 e[2][0]=0;e[2][1]=1;
 e[3][0]=-1;e[3][1]=0;
 e[4][0]=0;e[4][1]=-1;
 e[5][0]=1;e[5][1]=1;
 e[6][0]=-1;e[6][1]=1;
 e[7][0]=-1;e[7][1]=-1;
 e[8][0]=1;e[8][1]=-1;

 for(i=0;i<IMAX;i++)
  for(j=0;j<IMAX;j++)
   for(k=0;k<9;k++)
    {
     temp=(e[k][0]*u[i][j])+(e[k][1]*v[i][j]);     
     feq[k][i][j]=w[k]*rho[i][j]*(1+(3.0*temp)+(4.5*temp*temp)-(1.5*ur[i][j]));
    }

  for(i=0;i<IMAX;i++)
  for(j=0;j<IMAX;j++)
   for(k=0;k<9;k++)
   f[k][i][j]=feq[k][i][j];
}


iteration()
{
 int i,j,k;
 t=t+dt;

 /*------------------collision-------------------*/

 for(i=0;i<IMAX;i++)
  for(j=0;j<(IMAX-1);j++)
   for(k=0;k<9;k++)
    f[k][i][j]=f[k][i][j]*(1-(1/tou))+(feq[k][i][j]/tou);
    
 /*---------------streaming--------------------*/
 for(i=0;i<IMAX;i++)
  for(j=0;j<(IMAX-1);j++)
    f[0][i][j]=f[0][i][j];

 for(i=IMAX-1;i>0;i--)
  for(j=0;j<(IMAX-1);j++)
   f[1][i][j]=f[1][i-1][j];

 for(i=0;i<IMAX;i++)
  for(j=IMAX-2;j>0;j--)
   f[2][i][j]=f[2][i][j-1];

 for(i=0;i<(IMAX-1);i++)
  for(j=0;j<(IMAX-1);j++)
   f[3][i][j]=f[3][i+1][j];

 for(i=0;i<IMAX;i++)
  for(j=0;j<(IMAX-1);j++)
   f[4][i][j]=f[4][i][j+1];

 for(i=IMAX-1;i>0;i--)
  for(j=IMAX-2;j>0;j--)
   f[5][i][j]=f[5][i-1][j-1];

 for(i=0;i<(IMAX-1);i++)
  for(j=IMAX-2;j>0;j--)
   f[6][i][j]=f[6][i+1][j-1];
 
 for(i=0;i<(IMAX-1);i++)
  for(j=0;j<(IMAX-1);j++)
   f[7][i][j]=f[7][i+1][j+1];

 for(i=IMAX-1;i>0;i--)
  for(j=0;j<(IMAX-1);j++)
   f[8][i][j]=f[8][i-1][j+1];



/*----------Boundary---------------*/

 i=0;
 for(j=1;j<(IMAX-1);j++)
 {
  f[1][i][j]=f[3][i][j];
  f[5][i][j]=((1-sig)*f[6][i][j])+(sig*f[7][i][j]);
  f[8][i][j]=((1-sig)*f[7][i][j])+(sig*f[6][i][j]);
 }

 i=IMAX-1;
 for(j=1;j<(IMAX-1);j++)
 {
  f[3][i][j]=f[1][i][j];
  f[7][i][j]=((1-sig)*f[8][i][j])+(sig*f[5][i][j]);
  f[6][i][j]=((1-sig)*f[5][i][j])+(sig*f[8][i][j]);
 }

 j=0;
 for(i=1;i<(IMAX-1);i++)
 {
  f[5][i][j]=((1-sig)*f[8][i][j])+(sig*f[7][i][j]);
  f[2][i][j]=f[4][i][j];
  f[6][i][j]=((1-sig)*f[7][i][j])+(sig*f[8][i][j]);
 }

 j=0;i=0;
 {
  f[1][i][j]=f[3][i][j];
  f[5][i][j]=f[7][i][j];
  f[2][i][j]=f[4][i][j];
  f[6][i][j]=feq[6][i][j];
  f[8][i][j]=feq[8][i][j];

 }

 j=0;i=IMAX-1;
 {
  f[3][i][j]=f[1][i][j];
  f[6][i][j]=f[8][i][j];
  f[2][i][j]=f[4][i][j];
  f[5][i][j]=feq[5][i][j];
  f[7][i][j]=feq[7][i][j];
 }
 
 j=IMAX-1;
 for(i=0;i<IMAX;i++)
  for(k=0;k<9;k++)
   f[k][i][j]=feq[k][i][j];
}

variable()
{

 int i,j,k;
 float temp;
 
 for(i=0;i<IMAX;i++)
 {
  for(j=0;j<(IMAX-1);j++)
   {
    rho[i][j]=un[i][j]=vn[i][j]=0.0;
    
    for(k=0;k<9;k++)
     rho[i][j] = rho[i][j] + f[k][i][j];
 
    for(k=0;k<9;k++)
    {
     un[i][j] = un[i][j] + (f[k][i][j]*e[k][0]);
     vn[i][j] = vn[i][j] + (f[k][i][j]*e[k][1]);
    }
     un[i][j]=un[i][j]/rho[i][j];
     vn[i][j]=vn[i][j]/rho[i][j];
     
     ur[i][j]=(un[i][j]*un[i][j])+(vn[i][j]*vn[i][j]);
   }
  rho[i][j]=1.0;
  un[i][j]=U;
  vn[i][j]=0;
  ur[i][j]=U*U;
 }


 for(i=0;i<IMAX;i++)
  for(j=0;j<(IMAX-1);j++)
   for(k=0;k<9;k++)

    {
     temp=(e[k][0]*un[i][j])+(e[k][1]*vn[i][j]);
     feq[k][i][j]=w[k]*rho[i][j]*(1.0+(3.0*temp)+(4.5*temp*temp)-(1.5*ur[i][j]));
    }
}

error()
{
 int i,j;
 float err1,err2; 
 err=err1=err2=0;
  its++;
 for(i=0;i<IMAX;i++)
  for(j=0;j<(IMAX-1);j++)
      {
	  err1 = err1 + (pow((un[i][j]-u[i][j]),2)+pow((vn[i][j]-v[i][j]),2));
          err2 = err2 + ur[i][j];
          
      }
 err=sqrt(err1)/sqrt(err2);
  for(i=0;i<IMAX;i++)
  for(j=0;j<(IMAX-1);j++)
  {
   u[i][j]=un[i][j];
   v[i][j]=vn[i][j];
  }
if(its%100==0)
printf(" err = %.12lf    %d\n",err,its);
}

output()
{
 int i,j;
 FILE *fo,*fp,*fs;

 fp= fopen("(VvsX)Kn.0.1&sig0.5.plt","w");

// fprintf(fp,"\n\t\tTITLE=\"2D\"\n\tVARIABLES=\"X\",\"Y\",\"U\",\"V\"\n\t ZONE T=\"BLOCK1\",I=%d,J=%d,F=POINT\n\n",IMAX,IMAX);
 fprintf(fp," x      y   \n");
  for(i=0;i<IMAX;i++)
//   for(j=0;j<(IMAX-1);j++)
      fprintf(fp,"\n  %f \t %.24lf  ",1.0*i/IMAX,v[i][150]/U);
  
fclose(fp);

fs= fopen("(vortex)Kn.0.1&sig0.5.plt","w");

 fprintf(fs,"\n\t\tTITLE=\"2D\"\n\tVARIABLES=\"X\",\"Y\",\"U\",\"V\"\n\t ZONE T=\"BLOCK1\",I=%d,J=%d,F=POINT\n\n",IMAX-1,IMAX-1);
// fprintf(fs," x      y   \n");
  for(i=0;i<IMAX;i++)
   for(j=0;j<(IMAX-1);j++)
      fprintf(fs,"\n  %f \t %f \t %.24lf \t %.24lf  ",1.0*i/IMAX,1.0*j/IMAX,u[i][j]/U,v[i][j]/U);

fclose(fs);

fo= fopen("(UvsY)Kn.0.1&sig0.5.plt","w");

// fprintf(fo,"\n\t\tTITLE=\"2D\"\n\tVARIABLES=\"X\",\"Y\",\"U\",\"V\"\n\t ZONE T=\"BLOCK1\",I=%d,J=%d,F=POINT\n\n",IMAX,IMAX);
 fprintf(fo," x      y   \n");
  //for(i=0;i<IMAX;i++)
   for(j=0;j<(IMAX-1);j++)
      fprintf(fo,"\n  %.24lf \t %f  ",u[150][j]/U,1.0*j/IMAX);

fclose(fo);
} 


main()
{
 its=0;
 initialize();
 do
 {
  iteration();
  variable();
  error();
  if(its % 500 ==0)
  output();
  }while(its < 18000 );
output();
}

