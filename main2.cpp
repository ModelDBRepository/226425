#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "time.h"

const int N=9;
const double dt=0.001;
const double E_na=50.,E_k=-100.,E_l=-67.,E_i=-80.,E_ca=+120,V_shp = 2.5, V_lth=-25.;
//const double g_ahp=2.5,g_na=100.,g_k=80.,g_l=0.2,g_Ca=1.0,a_ca = 0.000025,t_ca=2000.,a_s=5.,b_s=0.5;
const double g_ahp=0.4,g_na=100.,g_k=80.,g_l=0.2,g_Ca=1.0,a_ca = 0.0002,t_ca=1000.,a_s=10.,b_s=0.5;
const double k_d=1.0;
int N_step1=2000000;
int N_step2=4000000;
double V[N],m[N],h[N],n[N],Ca[N],s[N];
double k_v[4][N],k_h[4][N],k_n[4][N],k_m[4][N],k_Ca[4][N],k_s[4][N];
double I_ext[N],divisor,tim=0;
double G_i[N*N];
double G_e[N*N];
FILE *dat1 = fopen("data1_1","wb");
FILE *dat2 = fopen("data2_1","wb");
FILE *dat3 = fopen("data3_1","wb");
double noise_int=0.2;
double sens_resp(double x)
{
	return 3.0/(1.+exp(-1.0*(x-6.)));
}
double gauss_gen()
{
	double rsq,x1,x2;
	do {
		x1 = -1.+2.*((double)rand())/((double)RAND_MAX);
		x2 = -1.+2.*((double)rand())/((double)RAND_MAX);
		rsq = x1*x1+x2*x2;
	}while((rsq>=1.0) || (rsq == 0.0));
	return x1*sqrt(-2.*log(rsq)/rsq);

}
double S(double x)
{
	return 1./(1+exp(-100*(x-20)));
}
double F_v(int num,double* V, double *k_v,double *s, double *k_s,double m,double h, double n, double Ca,double I_ext)
{
	double I_k = g_k*n*n*n*n*(V[num]+k_v[num]*divisor-E_k);
	double I_k_Ca = g_ahp*Ca/(Ca+k_d)*(V[num]+k_v[num]*divisor-E_k);
	double I_Ca = g_Ca/ ( 1.+exp( -(V[num]+k_v[num]*divisor - V_lth)/V_shp ) ) * (V[num]+k_v[num]*divisor-E_ca);
	double I_na = g_na*m*m*m*h*(V[num]+k_v[num]*divisor-E_na);
	double I_leak = g_l*(V[num]+k_v[num]*divisor-E_l);
	double I_syn=0.0;
	for (int i=0;i<N;i++)
	{
		I_syn=I_syn+G_i[num*N+i]*(s[i]+k_s[i]*divisor)*(-E_i+V[num]+k_v[num]*divisor);
	}
	double res = -(I_k+I_k_Ca+I_Ca+I_na+I_leak+I_syn+I_ext);
	return res;
}
double F_m(double V,double m)
{
	double a_m = 0.32*(54.+V)/ ( -exp(-(V+54.)*0.25) +1. ); 
	double b_m = 0.28*(V+27.)/ ( exp((V+27)*0.2) -1. );
	return a_m*(1-m)-b_m*m;

}
double F_h(double V,double h)
{
	double a_h = 0.128*( exp(-(50.+V)/18) ); 
	double b_h = 4./ ( exp(-(V+27)*0.2) +1. );
	return a_h*(1-h)-b_h*h;

}
double F_n(double V,double n)
{
	double a_n = 0.032*(52.+V)/ ( -exp(-(52.+V)*0.2) +1. ); 
	double b_n = 0.5* exp((-57.-V)/40. );
	return a_n*(1-n)-b_n*n;

}
double F_s(double V,double s)
{
	return a_s*S(V)*(1-s)-b_s*s;

}
double F_Ca(double V,double Ca)
{
	return -a_ca*g_Ca/ ( 1.+exp( -(V - V_lth)/V_shp ) ) * (V-E_ca) - Ca/t_ca;	
}
void solve(int Tau, int flag)
{
	
	for (int counter=0;counter<Tau;counter++)
	{
		tim+=dt;
		if (counter%100000==0) printf("%g \n",(double)counter/(double)Tau);
		if ((counter%10==0)&&(flag))
		{
			fwrite(V,sizeof(double)*N,1,dat1);
			fwrite(s,sizeof(double)*N,1,dat3);
		}
	/*	if (counter==440000)
		{
			fprintf(dat2,"V[1]=%g;m[1]=%g;h[1]=%g;n[1]=%g;m_T[1]=%g;h_T[1]=%g;Ca[1]=%g; \n",V[1],m[1],h[1],n[1],m_T[1],h_T[1],Ca[1]);
			fprintf(dat2,"V[2]=%g;m[2]=%g;h[2]=%g;n[2]=%g;m_T[2]=%g;h_T[2]=%g;Ca[2]=%g; \n",V[2],m[2],h[2],n[2],m_T[2],h_T[2],Ca[2]);
		}
		if (counter==470000)
		{
			fprintf(dat2,"V[0]=%g;m[0]=%g;h[0]=%g;n[0]=%g;m_T[0]=%g;h_T[0]=%g;Ca[0]=%g; \n",V[0],m[0],h[0],n[0],m_T[0],h_T[0],Ca[0]);
			fprintf(dat2,"V[3]=%g;m[3]=%g;h[3]=%g;n[3]=%g;m_T[3]=%g;h_T[3]=%g;Ca[3]=%g; \n",V[3],m[3],h[3],n[3],m_T[3],h_T[3],Ca[3]);
		}*/
		divisor=0.;
		for (int i=0;i<N;i++)
		{
			k_v[0][i] = F_v(i,V,k_v[0],s,k_s[0],m[i],h[i],n[i],Ca[i],I_ext[i])*dt;
			k_m[0][i] = F_m(V[i],m[i])*dt;
			k_h[0][i] = F_h(V[i],h[i])*dt;
			k_n[0][i] = F_n(V[i],n[i])*dt;
			k_s[0][i] = F_s(V[i],s[i])*dt;
			k_Ca[0][i] = F_Ca(V[i],Ca[i])*dt;
		}
			
		divisor=0.5;	
		for (int i=0;i<N;i++)
		{
			k_v[1][i] = F_v(i,V,k_v[0],s,k_s[0],m[i]+k_m[0][i]/2.0,h[i]+k_h[0][i]/2.0,n[i]+k_n[0][i]/2.0,Ca[i]+k_Ca[0][i]/2.0,I_ext[i])*dt;
			k_m[1][i] = F_m(V[i]+k_v[0][i]/2.0,m[i]+k_m[0][i]/2.0)*dt;
			k_h[1][i] = F_h(V[i]+k_v[0][i]/2.0,h[i]+k_h[0][i]/2.0)*dt;
			k_n[1][i] = F_n(V[i]+k_v[0][i]/2.0,n[i]+k_n[0][i]/2.0)*dt;
			k_s[1][i] = F_s(V[i]+k_v[0][i]/2.0,s[i]+k_s[0][i]/2.0)*dt;
			k_Ca[1][i] = F_Ca(V[i]+k_v[0][i]/2.0,Ca[i]+k_Ca[0][i]/2.)*dt;
		}
			
			
		for (int i=0;i<N;i++)
		{
			k_v[2][i] = F_v(i,V,k_v[1],s,k_s[1],m[i]+k_m[1][i]/2.0,h[i]+k_h[1][i]/2.0,n[i]+k_n[1][i]/2.0,Ca[i]+k_Ca[1][i]/2.0,I_ext[i])*dt;
			k_m[2][i] = F_m(V[i]+k_v[1][i]/2.0,m[i]+k_m[1][i]/2.0)*dt;
			k_h[2][i] = F_h(V[i]+k_v[1][i]/2.0,h[i]+k_h[1][i]/2.0)*dt;
			k_n[2][i] = F_n(V[i]+k_v[1][i]/2.0,n[i]+k_n[1][i]/2.0)*dt;
			k_s[2][i] = F_s(V[i]+k_v[1][i]/2.0,s[i]+k_s[1][i]/2.0)*dt;
			k_Ca[2][i] = F_Ca(V[i]+k_v[1][i]/2.0,Ca[i]+k_Ca[1][i]/2.)*dt;
		}

		divisor=1.;
		for (int i=0;i<N;i++)
		{
			k_v[3][i] = F_v(i,V,k_v[2],s,k_s[2],m[i]+k_m[2][i],h[i]+k_h[2][i],n[i]+k_n[2][i],Ca[i]+k_Ca[2][i],I_ext[i])*dt;
			k_m[3][i] = F_m(V[i]+k_v[2][i],m[i]+k_m[2][i])*dt;
			k_h[3][i] = F_h(V[i]+k_v[2][i],h[i]+k_h[2][i])*dt;
			k_n[3][i] = F_n(V[i]+k_v[2][i],n[i]+k_n[2][i])*dt;
			k_s[3][i] = F_s(V[i]+k_v[2][i],s[i]+k_s[2][i])*dt;
			k_Ca[3][i] = F_Ca(V[i]+k_v[2][i],Ca[i]+k_Ca[2][i]/2.)*dt;
		}
//--------------here edited

		for (int i=0;i<N;i++)
		{
			V[i] += (k_v[0][i]+2.0*k_v[1][i]+2.0*k_v[2][i]+k_v[3][i])/6.0;
			m[i] += (k_m[0][i]+2.0*k_m[1][i]+2.0*k_m[2][i]+k_m[3][i])/6.0;
			h[i] += (k_h[0][i]+2.0*k_h[1][i]+2.0*k_h[2][i]+k_h[3][i])/6.0;
			n[i] += (k_n[0][i]+2.0*k_n[1][i]+2.0*k_n[2][i]+k_n[3][i])/6.0;
			s[i] += (k_s[0][i]+2.0*k_s[1][i]+2.0*k_s[2][i]+k_s[3][i])/6.0;
			Ca[i] += (k_Ca[0][i]+2.0*k_Ca[1][i]+2.0*k_Ca[2][i]+k_Ca[3][i])/6.0;

			V[i] = V[i] + noise_int*sqrt(dt)*gauss_gen();
		}
	}
}
int main()
{
	srand(time(NULL));
	double Delta=2.0;
	double sens_stim[N];
	
	
	for (int i=0;i<N;i++)
	{
		m[i]=n[i]=h[i]=Ca[i]=0.01;
		Ca[i]=0.0001;
		V[i] = -60.6+0.0*rand()/RAND_MAX;
		I_ext[i]=-0.0-0.*(double)rand()/(double)RAND_MAX;
	}
	for(int i=0;i<N;i++)
		for (int j=0;j<N;j++) 
		{
			G_i[i*N+j]=20.0;
			G_e[i*N+j]=0.0;
		}
	for(int i=0;i<N;i++)
		 G_i[i*N+i]=0.0;

	solve(N_step1,1);
	for (int i=0;i<N;i++)
	{
		double shift = -2.+4.*(double)i/(double)(N-1);
		//sens_stim[i] = 4.0+(Delta*(double)(N-1.-i)/(double)(N-1));
		sens_stim[i] = 2.0;
		I_ext[i]=-0.0-sens_resp(sens_stim[i]-shift);
		//I_ext[i]=-sens_resp(sens_stim[i]-shift);
	}
	fwrite(I_ext,sizeof(double)*N,1,dat2);
	fwrite(sens_stim,sizeof(double)*N,1,dat2);

	solve(N_step2,1);




	return 0;
}
