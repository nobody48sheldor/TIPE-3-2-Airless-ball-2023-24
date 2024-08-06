#include <iostream>
#include <format>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

using namespace std;


// fonctions generales

auto linspace(double lower_bound, double upper_bound, int taille)
{
	vector<double> array;
	array.resize(taille);

	for (int i=0; i<taille; i++)
	{
		array[i] = lower_bound + (i*(upper_bound - lower_bound)/(taille-1));
	}
	return array;
}

double partial_derivative_1(double u_, double u, double differential)
{
	return ((u_ - u)/differential);
}

double partial_derivative_2(double u__, double u_, double u, double differential)
{
	return ((u__ - 2*u_ + u)/(differential*differential));
}

// fonction d ecriture des donnees


int write_data(int n, int nt, vector<vector<vector<double>>> u_x, vector<vector<vector<double>>> u_y)
{
	ofstream output_file;
	output_file.open( "data/data_x.txt" );
	for (int t=0; t<nt; t++)
	{
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				output_file << u_x[t][i][j] << '/';
			}
		}
	}
	output_file.close();

	output_file.open( "data/data_x.txt" );
	for (int t=0; t<nt; t++)
	{
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				output_file << u_y[t][i][j] << '/';
			}
		}
	}

	return 0;
}

// fonction d'initialisation

auto initialisation(int n, int nt, vector<double> X, vector<double> Y, double v, double R, double g, double rho)
{
	vector<vector<double>> u_x_init;
	vector<vector<double>> u_y_init;

	// on alloue la memoire pour nos vecteurs

	u_x_init.resize(n);
	u_y_init.resize(n);
	for (int i=0; i<n; i++)
	{
		u_x_init[i].resize(n);
		u_x_init[i].resize(n);
	}


	// on donne la foncrion d'intialisation pour le flux AUTOUR d'une balle avec un nombre de Reynold << 1

	double r;
	double theta;
	double u_r;
	double u_theta;
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
		{
			theta = atan(Y[j]/X[i]);
			r = sqrt( X[i]*X[i] + Y[j]*Y[j] );

			if ( (r*r) < (R*R) )
			{
				u_r = -v * cos(theta) * (1 - (3*R)/(2*r) + (R*R*R/(2*r*r*r)));
				u_theta = -v * sin(theta) * (1 - (3*R)/(4*r) - (R*R*R/(4*r*r*r)));
				u_x_init[i][j] = cos(theta) * u_r + sin(theta) * u_theta;
				u_y_init[i][j] = sin(theta) * u_r - cos(theta) * u_theta;
			}
			else
			{	
				u_x_init[i][j] = 0;
				u_x_init[i][j] = 0;
			}
		}
	}
	return(u_x_init, u_y_init);
}



// fonction de résolution

	// 
auto solve(int n, int nt, int t, double dx, double dy, double dt, vector<vector<vector<double>>> u_x, vector<vector<vector<double>>> u_y, vector<double> X, vector<double> Y, double v, double R, double g, double rho)
{
	vector<vector<double>> u_x_new;
	vector<vector<double>> u_y_new;
	u_x_new.resize(n);
	u_y_new.resize(n);

	// on alloue la memoire pour nos vecteurs
	
	u_x_new.resize(n);
	u_y_new.resize(n);
	for (int i=0; i<n; i++)
	{
		u_x_new[i].resize(n);
		u_x_new[i].resize(n);
	}

	// on resout sans oublier les conditions aux bords

	for (int j=0; j<n; j++)
	{
		u_x[t+1][0][j] = v;
		u_y[t+1][0][j] = 0;
	}
	for (int i=1; i<(n-1); i++)
	{
		u_x[t+1][i][0] = v;
		u_y[t+1][i][n-1] = 0;
		for (int j=1; j<(n-1); j++)
		{
			double r = sqrt( X[i]*X[i] + Y[j]*Y[j] );
			if ( (r*r)< (R*R) )
			{
				u_x[t+1][i][j] = u_x[t][i][j] + (
						-0.5* partial_derivative_1(u_x[t][i+1][j]*u_x[t][i+1][j] + u_y[t][i+1][j]*u_y[t][i+1][j], u_x[t][i][j]*u_x[t][i][j] + u_y[t][i][j]*u_y[t][i][j], dx)
						- u_x[t][i][j] * partial_derivative_1(u_x[t][i+1][j], u_x[t][i][j], dx)
						- u_y[t][i][j] * partial_derivative_1(u_x[t][i][j+1], u_x[t][i][j], dy)
						) * dt;




				u_y[t+1][i][j] = u_y[t][i][j] + (
						-0.5* partial_derivative_1(u_x[t][i][j+1]*u_x[t][i][j+1] + u_y[t][i][j+1]*u_y[t][i][j+1], u_x[t][i][j]*u_x[t][i][j] + u_y[t][i][j]*u_y[t][i][j], dy)
						- u_x[t][i][j] * partial_derivative_1(u_y[t][i+1][j], u_y[t][i][j], dx)
						- u_y[t][i][j] * partial_derivative_1(u_y[t][i][j+1], u_y[t][i][j], dy)
						-g
						) * dt;


			}
			else
			{
				u_x[t+1][i][j] = 0;
				u_y[t+1][i][j] = 0;
			}

		}
	}
	return 0;
}



// fonctions main

int main()
{
	// definitions des valeurs
	
	int n = 100;
	int nt = 100;

	double g = 9.81;
	double rho = 1.204;
	double R = 0.2;
	double v = 1;
	vector<double> X = linspace(-0.5, 0.5, n);
	vector<double> Y = linspace(-0.5, 0.5, n);
	vector<double> T = linspace(0.0, 3.0, n);

	vector<vector<vector<double>>> u_x;
	vector<vector<vector<double>>> u_y;
	u_x.resize(nt);
	u_y.resize(nt);


	// on alloue la memoire pour nos vecteurs
	
	for (int t=0; t<nt; t++)
	{
		u_x[t].resize(n);
		u_y[t].resize(n);
		for (int i=0; i<n; i++)
		{
			u_x[t][i].resize(n);
			u_y[t][i].resize(n);
		}
	} // on a u_x et u_y de taille nt * n * n
	
	auto vect = initialisation(n, nt, X, Y, v, R, g, rho);

	u_x[0], u_y[0] = vect;

	double dx = X[1] - X[0];
	double dy = Y[1] - Y[0];
	double dt = T[1] - T[0];


	
	return 0;

	//for (int t=0; t<(nt-1); t++)
	//{
		//cout << t << endl;
		//solve(n, nt, t, dx, dy, dt, u_x, u_y, X, Y, v, R, g, rho);
	//}

	//write_data(n, nt, u_x, u_y);
	//return 0;
}
