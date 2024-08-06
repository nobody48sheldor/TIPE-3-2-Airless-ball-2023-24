#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

using namespace std;


// general function

auto linspace(double low, double up, int n)
{
	vector<double> array;
	array.resize(n);

	for (int i=0; i<n; i++)
	{
		array[i] = low + (i*(up-low)/(n-1));
	}
	return array;
}

double grad(double u__,double u_,double u,double dx)
{
	return( ((u_ - u__) + (u - u_))/(2*dx) );
}

double grad2(double u__,double u_,double u,double dx)
{
	return( (u - 2*u_ + u__ )/(dx*dx)  );
}

int write_data(vector<double> data, int n)
{
  ofstream output_file;
  output_file.open("data/data_U_lin.txt");
  for (int i=0; i < n; i++)
  {
    output_file << data[i] << '/';
  }
  output_file.close();
  return 0;
} 


// init function

auto init_func(double x, double y)
{
	double v_x = x+y;
	double v_y = x-y;
	vector<double> v;
	v.resize(2);
	v[0] = v_x;
	v[1] = v_y;
	return(v);
}


// solve

auto solve(vector<vector<vector<vector<double>>>> u, vector<vector<vector<double>>> U, vector<vector<vector<double>>> P, double g, double mu,double rho, double constante, vector<double> y, double dx, double dy, double dt, int n, int nt)
{
	for(int time=0;time<nt;time++)
	{
		cout << "time = " << time << endl;
		u[time+1][0] = u[0][0];
		u[time+1][n] = u[0][n];

		// border
		for (int xpos=0; xpos<(n+1); xpos++)
		{
			U[time+1][0][xpos] = sqrt( u[0][0][xpos][0]*u[0][0][0][0] + u[time+1][0][xpos][1]*u[0][0][xpos][1] );
			U[time+1][n][xpos] = sqrt( u[0][n][xpos][0]*u[0][n][0][0] + u[time+1][0][xpos][1]*u[0][n][xpos][1] );
		}
		// border
		for(int ypos=1;ypos<(n);ypos++)
		{
			cout << "ypos = " << ypos << endl;
			u[time+1][ypos][0][0] = u[0][ypos][0][0];
			u[time+1][ypos][0][1] = u[0][ypos][0][0];
			u[time+1][ypos][n][0] = u[0][ypos][0][1];
			u[time+1][ypos][n][1] = u[0][ypos][0][1];
			U[time+1][ypos][0] = sqrt( u[0][ypos][0][0]*u[0][ypos][0][0] + u[time+1][ypos][0][1]*u[0][ypos][0][1] );
			U[time+1][ypos][n] = sqrt( u[0][ypos][0][0]*u[0][ypos][0][0] + u[time+1][ypos][0][1]*u[0][ypos][0][1] );
			for(int xpos=1;xpos<(n);xpos++)
			{
				u[time+1][ypos][xpos][0] = u[time][ypos][xpos][0] + ( -1*u[time][ypos][xpos][0] * ( grad(U[time][ypos][xpos+1],U[time][ypos][xpos], U[time][ypos][xpos-1], dx) + grad(U[time][ypos+1][xpos], U[time][ypos][xpos], U[time][ypos-1][xpos], dy) ) + mu * ( grad2(U[time][ypos][xpos+1], U[time][ypos][xpos], U[time][ypos][xpos-1], dx) + grad2(U[time][ypos+1][xpos], U[time][ypos][xpos], U[time][ypos-1][xpos], dy) ) - (1.0/rho) * ( grad(P[time][ypos][xpos+1], P[time][ypos][xpos], P[time][ypos][xpos-1], dx) + grad(P[time][ypos+1][xpos], P[time][ypos-1][xpos], P[time][ypos-1][xpos], dx) ) )*dt;

				u[time+1][ypos][xpos][1] = u[time][ypos][xpos][1] + ( -1*u[time][ypos][xpos][0] * ( grad(U[time][ypos][xpos+1],U[time][ypos][xpos], U[time][ypos][xpos-1], dx) + grad(U[time][ypos+1][xpos], U[time][ypos][xpos], U[time][ypos-1][xpos], dy) ) + mu * ( grad2(U[time][ypos][xpos+1], U[time][ypos][xpos], U[time][ypos][xpos-1], dx) + grad2(U[time][ypos+1][xpos], U[time][ypos][xpos], U[time][ypos-1][xpos], dy) ) - (1.0/rho) * ( grad(P[time][ypos][xpos+1], P[time][ypos][xpos], P[time][ypos][xpos-1], dx) + grad(P[time][ypos+1][xpos], P[time][ypos-1][xpos], P[time][ypos-1][xpos], dx) ) - g)*dt;

				U[time+1][ypos][xpos] = sqrt( u[time+1][ypos][xpos][0]*u[time+1][ypos][xpos][0] + u[time+1][ypos][xpos][1]*u[time+1][ypos][xpos][1] );

				P[time+1][ypos][xpos] = constante - 0.5*rho*U[time+1][ypos][xpos]*U[time+1][ypos][xpos] - rho*g*y[ypos];
			}
		}
	}
	return u, U, P;
}


// main function

int main()
{
	// init variables
	
	int n = 100;
	int nt = 100;

	double g = 9.81;
	double mu = 9.81;
	double rho = 1.2;
	double v = 5;
	double P_atm = 101300;
	double constante = 0.5*rho*v*v + P_atm;
	vector<double> x = linspace(-1.0, 1.0, n);
	vector<double> y = linspace(-1.0, 1.0, n);
	vector<double> t = linspace(0.0, 10.0, n);

	vector<vector<vector<vector<double>>>> u;
	u.resize(n);
	for (int i=0; i<n; i++)
	{
		u[i].resize(n);
		for (int j=0; j<n; j++)
		{
			u[i][j].resize(n);
			for (int k=0; k<n; k++)
			{
				u[i][j][k].resize(n);
			}
		}
	} // empty speed vect

	vector<vector<vector<double>>> U;
	U.resize(n);
	for (int i=0; i<n; i++)
	{
		U[i].resize(n);
		for (int j=0; j<n; j++)
		{
			U[i][j].resize(n);
		}
	} // empty speed scalar field

	vector<vector<vector<double>>> P;
	P.resize(n);
	for (int i=0; i<n; i++)
	{
		P[i].resize(n);
		for (int j=0; j<n; j++)
		{
			P[i][j].resize(n);
		}
	} // empty Pressure scalar field


	double dx = x[1] - x[0];
	double dy = y[1] - y[0];
	double dt = t[1] - t[0];


	// solve
	
	// init
	for (int ypos=0; ypos<n; ypos++)
	{
		for (int xpos=0; xpos<n; xpos++)
		{
			u[0][ypos][xpos][0] = init_func(x[xpos],y[ypos])[0];
			u[0][ypos][xpos][1] = init_func(x[xpos],y[ypos])[1];
			U[0][ypos][xpos] = sqrt( (u[0][ypos][xpos][0]*u[0][ypos][xpos][0] + u[0][ypos][xpos][1]*u[0][ypos][xpos][1]) );
			P[0][ypos][xpos] = constante - 0.5*rho*U[0][ypos][xpos]*U[0][ypos][xpos] - rho*g*y[ypos];
		}
	}

	u, U, P = solve(u, U, P, g, mu, rho, constante, y, dx, dy, dt,n,nt);
	// writing U
	vector<double> U_lin;
	U_lin.resize(n*n*nt);
	for (int time=0; time<nt; time++)
	{
		for (int ypos=0; ypos<n; ypos++)
		{
			for (int xpos=0; xpos<n; xpos++)
			{
				U_lin[xpos + n*ypos + n*nt*time] = U[time][ypos][xpos];
			}
		}
	};

	write_data(U_lin, n*n*nt);
		
	return 0;
}
