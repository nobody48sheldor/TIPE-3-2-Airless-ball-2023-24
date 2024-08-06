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

int write_data(vector<double> data, int n, int index)
{
  ofstream output_file;
  output_file.open("data/data_" + to_string(index) + ".txt");
  for (int i=0; i < n; i++)
  {
    output_file << data[i] << '/';
  }
  output_file.close();
  return 0;
}

int infos(double dt, double m, double r, double K_D, double K_M, double k_r, double J, double alpha, double v_0, vector<double> v_x , vector<double> v_z, double w_0)
{
	cout << "dt = " << dt << endl;
	cout << "masse = " << m << endl;
	cout << "rayon = " << r << endl;
	cout << "K_D = " << K_D << endl;
	cout << "K_M = " << K_M << endl;
	cout << "k_r = " << k_r << endl;
	cout << "J = " << J << endl;
	cout << "alpha = " << alpha*(180/3.141591628) << endl;
	cout << "v_0 = " << v_0 << endl;
	cout << "v_x[0] = " << v_x[0] << endl;
	cout << "v_z[0] = " << v_z[0] << endl;
	cout << "w_0 = " << w_0 << endl;
	return 0;
}



// solve

double solve_theta(double theta_time, double v_x_time_1, double v_x_time, double v_z_time_1, double v_z_time)
{
	// return ( theta_time + asin( ((v_x_time*v_z_time_1) - (v_x_time_1*v_z_time)) / ( sqrt( (v_x_time*v_x_time + v_z_time*v_z_time) *  ( v_x_time_1*v_x_time_1 + v_z_time_1*v_z_time_1 ) ) ) ) );
	return ( atan(v_z_time/v_x_time) );
}

double solve_speed_x(double v_x, double v_z, double t, double dt, double theta, double m, double g, double K_D, double K_M, double k_r, double J, double w_0)
{
	// cout << (-cos(theta) * K_D * ( v_x*v_x + v_z*v_z ) ) << endl;
	// cout << -( sin(theta) * K_M * (w_0*exp(-(k_r/J)*t)) * ( sqrt( v_x*v_x + v_z*v_z ) ) ) << endl;
	// cout << K_D << endl;
	// cout << K_M << endl;
	return( v_x + ( (-cos(theta) * K_D * ( v_x*v_x + v_z*v_z ) ) - ( sin(theta) * K_M * (w_0*exp(-(k_r/J)*t)) * ( sqrt( v_x*v_x + v_z*v_z ) ) ) )*(dt/m) );
}

double solve_speed_z(double v_x, double v_z, double t, double dt, double theta, double m, double g, double K_D, double K_M, double k_r, double J, double w_0)
{
	// cout << ( - m*g - ( sin(theta) * K_D * ( v_x*v_x + v_z*v_z ) ) + ( cos(theta) * K_M * ( w_0*exp(-(k_r/J)*t) ) *  ( sqrt(v_x*v_x + v_z*v_z) ) ) )*(dt/m)<< endl;
	return( v_z + ( - m*g - ( sin(theta) * K_D * ( v_x*v_x + v_z*v_z ) ) + ( cos(theta) * K_M * ( w_0*exp(-(k_r/J)*t) ) *  ( sqrt(v_x*v_x + v_z*v_z) ) ) )*(dt/m) );
}
	
// main function

int main(int argc, char *argv[])
{
	// init variables
	
	int nt = 50000;
	vector<double> t = linspace(0.0, 10.0, nt);
	double dt = t[1] - t[0];

	double g = 9.81;

	double m = 0.6;
	double r = 0.12;
	double K_D = 0.5 * 1.2 * 3.141591628 * r *r * 0.4;
	//double K_D = 0.5 * 1.2 * 3.141591628 * r *r * 0.54;
	// cout << K_D << endl;
	double K_M = 0.5 * 1.2 * 3.141591628 * r * r * r;
	double k_r = 0.001;
	double J = 0.667*m*r*r;

	vector<double> v_x;
	vector<double> v_z;
	vector<double> theta;
	v_x.resize(nt);
	v_z.resize(nt);
	theta.resize(nt);

	vector<double> x;
	vector<double> z;
	x.resize(nt);
	z.resize(nt);
	x[0] = 0;
	z[0] = 1.9;
	
	// cout << "ARGC = " << argc << endl;

	double alpha = 75*(3.141591628)/180;
	double v_0 = 9.5;


	if (argc == 3)
	{
		v_0 = strtod(argv[1], NULL)*1.0;
		alpha = strtod(argv[2], NULL)*(3.141591628)/180;
		// cout << v_0 << " v_0 C++" << endl;
		// cout << alpha*(180/3.141591628) << " alpha C++" << endl;
	};

	//cout << v_0 << endl;
	//cout << alpha << endl;


	double w_0 = 2*3.141591628*3;


	theta[0] = alpha;
	v_x[0] = cos(theta[0])*v_0;
	v_z[0] = sin(theta[0])*v_0;


	// infos(dt, m, r, K_D, K_M, k_r, J, alpha, v_0, v_x ,v_z, w_0);


	// solve

	for (int time=0; time<nt; time++)
	{
		v_x[time+1] = solve_speed_x(v_x[time], v_z[time], t[time], dt, theta[time], m, g, K_D, K_M, k_r, J, w_0);
		v_z[time+1] = solve_speed_z(v_x[time], v_z[time], t[time], dt, theta[time], m, g, K_D, K_M, k_r, J, w_0);
		theta[time+1] = solve_theta(theta[time], v_x[time+1], v_x[time], v_z[time+1], v_z[time]);

		//cout << "time = " << time << endl;
		//cout << "v_x = " << v_x[time] << endl;
		//cout << "v_z = " << v_z[time] << endl;
		//cout << "theta = " << (180/3.141591628)*theta [time] << endl;
		//cout << endl;

	}

	for (int time=0; time<nt; time++)
	{
		x[time+1] = x[time] + v_x[time]*dt;
		z[time+1] = z[time] + v_z[time]*dt;
		if (z[time+1] < 0)
		{
			z[time+1] = 0;
		}
	}

	write_data(x,nt,1);
	write_data(z,nt,2);

	return 0;
	return 0;
}
