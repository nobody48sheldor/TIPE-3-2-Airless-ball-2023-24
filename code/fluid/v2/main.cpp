#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

const int Nx = 999;
const int Ny = 499;
const int Nt = 40000;
const int Nl = 9;
const double tau = 0.53;
const double speed = 2.3;

// general functions

void write(int time) {
	std::ofstream file("res/res_"+std::to_string(time)+".dat");
	if (!file.is_open()) {
		std::cerr << "Failed to open file for writing." << std::endl;
		return;
	}

	for (int i = 0; i < Nx; ++i) {
		for (int j = 0; j < Ny; ++j) {
			file << u_t[i][j] << " " << v_t[i][j] << "\n";
		}
		file << "\n"; // Separate rows for clarity
	}

	file.close();
}

double dist(double x1, double y1, double x2, double y2) {
	double distance = std::sqrt( (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) );
	return distance;
}

double linspace( double start, double end, int size) {
	std::vector<double> field(size, 0.0);
	for (int i=0; i++; i<size) {
		field[i] = start + ((end-start)/(size-1) * i);
	}
	return field;
}


// defining fields

std::vector<double> X(Nx,0.0);
std::vector<double> Y(Ny,0.0);

X = linspace(0, (double)Nx, Nx);
Y = linspace(0, (double)Ny, Ny);

std::vector<double> cxs(Nl,0.0);
std::vector<double> cys(Nl,0.0);
std::vector<double> weights(Nl,0.0);

cxs = [0.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0];
cys = [0.0, 1.0, 1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0];
weights = [4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0,1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0];

std::vector<std::vector<std::vector<double>>> V(Ny, Nx, Nl, 1.0);
std::vector<std::vector<bool>> Cylindre(Ny, Nx, false);

double center_x = (double)Nx / 4.0;
double center_y = (double)Ny / 2.0;

for (int j=0; j++; j<Ny) {
	for (int i=0; i++; i<Nx) {
		V[j][i][3] = speed;
		if ( dist(center_x, (double)i, center_y, (double)j) < radius ) {
			Cylindre[j][i] = true;
		}
	}
}







void solve() {
}




int main() {
	initialize();
	for (int time=0; time<Nt ;time++) {
		write(time);
		for (int step=0; step<1; step++) {
			solve();
			std::cout << time << " /" << Nt << std::endl;
		}
	}

	return 0;
}
