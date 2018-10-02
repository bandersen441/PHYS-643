/* Bridget Andersen, 9/29/18
   PHYS 643, Computational Exercise #1
   This script contains all the code for simulating a white dwarf density profile using
   4th order Runge Kutta integration.
   Outputs:
   rk4_pc_####.txt : The density profiles for different values of central density p_c. 
   #### is a number indicating the order in which the files were created. A higher #### 
   corresponds to a higher p_c value.
   pc_mass_radius.txt : The mass-radius relation.
*/
#include <iostream>
#include <math.h>
#include <vector>
#include <fstream>
using namespace std;

/* Some useful constants */
double PI = M_PI;
double Y_e = 0.5;
double G = 6.67259 * pow(10., -8);
double h_bar = 1.05457266 * pow(10., -27);
double m_e = 9.1093897 * pow(10., -28);
double m_p = 1.6726231 * pow(10., -24);
double c = 2.99792458 * pow(10., 10);
double M_sun = 1.99 * pow(10., 33);
/* Formulas for calculating the constants K_nr and K_r in cgs units */
double K_nr = pow(h_bar, 2) * pow(3. * pow(PI, 2), 2. / 3.) * pow(Y_e, 5. / 3.) / (5. * m_e * pow(m_p, 5./3.));
double K_r = h_bar * c * pow(3. * pow(PI, 2), 1. / 3.) * pow(Y_e, 4. / 3.) / (4 * pow(m_p, 4./3.));

/* Calculates the equation of state for non-relativistic electrons */
double P_nr(double p) {
  return K_nr * pow(p, 5. / 3.);
}

/* Calculates the equation of state for relativistic electrons */
double P_r(double p) {
  return K_r * pow(p, 4. / 3.);
}

/* Calculate the Paczynksi (1983) equation of state */
double P(double p) {
  return pow(pow(P_nr(p), -2) + pow(P_r(p), -2), -0.5);
}

/* To simplify calculations, we define:
A = [5/3 * (P / P_nr)^2 + 4/3 * (P / P_r)^2]
This function calculates this quantity */
double A(double p) {
  return (5. / 3.) * pow(P(p) / P_nr(p), 2) + (4. / 3.) * pow(P(p) / P_r(p), 2);
}

/* Calculate the density derivative dp/dr given the density p, mass m, and radius r */
double dp_dr(double r, double p, double m) {
  return - G * m * pow(p, 2) / (P(p) * A(p) * pow(r, 2));
}

/* Calculate the mass derivative dm/dr give the density and radius r */
double dm_dr(double r, double p) {
  return 4 * PI * pow(r, 2) * p;
}

/* Calculate the Fermi Energy given a central density */
double E_F(double p) {
  double x = pow(p * Y_e / pow(10., 6), 1./3.);
  return m_e * pow(c, 2) * (sqrt(1 + pow(x, 2)) - 1);
}

/* Complete a 4th order Runge-Kutta integration of both the density and the mass as a
   function of radius, given an initial central density p_c. 
   Output in columns to the file: rk4_pc_(value of p_c).txt */
pair<double,double> rk4_integ(double p_c, bool save, int filenum) {
  /* Set up file output */
  ofstream fs;
  if(save) {
	char c[10];
	sprintf(c, "%04d", filenum);
	string s(c);
	string filename = "./density_profiles/rk4_pc_" + s + ".txt";
	fs.open(filename);
	fs.precision(20);
  }

  /* Set the radius step size to: 1 km = 10^5 cm */
  unsigned long long dr = pow(10., 5);
  /* Set the initial radius to 1 cm */
  unsigned long long r_0 = 1;
  /* Set the final radius to 10^10 cm */
  unsigned long long r_f = 10000000000;
  /* Calculate the number of iterations needed to reach the final radius 
	 (ignore the inital radius in this determination) */
  unsigned long long n_iter = r_f / dr;
  /* Calculate the initial mass */
  double m_0 = 4 * PI * p_c * pow(r_0, 3)/ 3.;

  /* Set the current density variable */
  double p = p_c;
  /* Set the current mass variable */
  double m = m_0;
  /* Set the currrent radius variable */
  unsigned long long r = r_0;

  /* Output the initial values */
  if(save) {
	fs << r << '\t' << p << '\t' << m / M_sun << endl;
  }

  unsigned long long progress = 0;
  //for(long i = 0; i < n_iter; i++) {
  while(p > 1.) {
	/* Calculate the RK4 increments for density (k_i) and mass (l_i) 
	   Note that these must be calculated in a particular order. */
	double k_0 = dr * dp_dr(r, p, m);
	double l_0 = dr * dm_dr(r, p);
	double k_1 = dr * dp_dr(r + 0.5 * dr, p + 0.5 * k_0, m + 0.5 * l_0);
	double l_1 = dr * dm_dr(r + 0.5 * dr, p + 0.5 * k_0);
	double k_2 = dr * dp_dr(r + 0.5 * dr, p + 0.5 * k_1, m + 0.5 * l_1);
	double l_2 = dr * dm_dr(r + 0.5 * dr, p + 0.5 * k_1);
	double k_3 = dr * dp_dr(r + dr, p + k_2, m + l_2);
	double l_3 = dr * dm_dr(r + dr, p + k_2);

	/* Calculate the next values of the density and the mass */
	p = p + (1. / 6.) * (k_0 + 2. * k_1 + 2 * k_2 + k_3);
	double m_check = m + (1. / 6.) * (l_0 + 2. * l_1 + 2 * l_2 + l_3);
	if(!isnan(m_check))
	  m = m_check;

	/* Iterate the radius */
	r = r + dr;

	/* Output new values */
	if(save) 
	  fs << r << '\t' << p << '\t' << m / M_sun << endl;
  }

  fs.close();
  pair<double,double> mr = make_pair(m, r);
  return mr;
}

/* Takes the maximum and minimum values and returns a vector of N values spaced 
   logarithmically within the given range. */
vector<double> log_space(double min_val, double max_val, long N) {
  double min_log = log10(min_val);
  double max_log = log10(max_val);
  double log_incr = (max_log - min_log) / N;

  double val = min_val;
  double log_val = min_log;
  vector<double> vals(N);
  vals[0] = min_val;
  for(int i = 1; i < N; i++) {
	log_val += log_incr;
	val = pow(log_val, 10);
	vals[i] = val;
  }
  return vals;
}

/* Takes the maximum and minimum values and returns a vector of N values spaced
   linearly withing the given range */
vector<double> lin_space(double min_val, double max_val, long N) {
  double incr = (max_val - min_val) / N;
  vector<double> vals(N);
  vals[0] = min_val;
  double val = min_val;
  for(int i = 1; i < N; i++) {
	val += incr;
	vals[i] = val;
  }
  return vals;
}

int main() {
  cout << "An RK4 Integration of the White Dwarf Mass-Radius Relation " << endl;
  /* Create an output file for masses and radii */
  ofstream fs;
  string filename = "pc_mass_radius.txt";
  fs.open(filename);
  fs.precision(20);

  /* Number of central densities to sample */
  int n_sample = 5000;
  
  /* Create a vector of central densities to iterate through */
  vector<double> p_c = log_space(pow(10,1), pow(10,20), n_sample);

  int progress = 0;
  for(int i = 0; i < p_c.size(); i++) {
	/* Start the integration */
	pair<double,double> mr;
	if(progress % 10 == 0)
	  mr = rk4_integ(p_c[i], 1, progress);
	else
	  mr = rk4_integ(p_c[i], 0, progress);
	double mass = mr.first;
	double radius = mr.second;
	double E_fermi = E_F(p_c[i]);
	/* Calculate the polytropic index */
	double ln_P = log(P(p_c[i]));
	double ln_p = log(p_c[i]);
	fs << p_c[i] << '\t' << mass / M_sun << '\t' << radius << '\t' << E_fermi << '\t' << ln_P << '\t' << ln_p << endl;

	progress += 1;
	if(progress % 100 == 0)
	  cout << "Progress: " << progress << " out of " << n_sample << " p_c's sampled..." << endl;
  }
  
  fs.close();
  return 0;
}
