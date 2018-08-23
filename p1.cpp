/*

s-wave scattering of massless scalar field by black hole 

compile with:

   g++ -std=c++11 -g -Wall -O2 p1.cpp -o p1 -lbbhutil

input parameters in terminal as:
(do NOT include .sdf extension in outfile)

 ./p1 <outfile> <lastpt> <save_pt> <nsteps> <save_step>
      <lam> <r2m> <rmin> <rmax> <dspn> <tol> <maxit>
      <ic_Dsq> <ic_r0> <ic_Amp> <check_step> <zero_pi> 
      <sommerfeld> <dspn_bound> <write_ires> <write_res>
      <write_itn> <hold_const>

where the value of any <parameter> can be set with the
following ordered pair of command line arguments:

 -parameter parameter_value

default values can be found at the start of main()
*/

#include <iostream>
#include <algorithm> // for max_element()
#include <vector> // for everything
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include <cmath> // for ICs
#include "bbhutil.h" // for output to .sdf
using namespace std;

#define _USE_MATH_DEFINES
#ifndef M_PI
cerr << "NO PI\n";
const double M_PI = 4.0*atan(1.0);
#endif

// **********************************************************
// **********************************************************
//                       FUNCTIONS
// **********************************************************
// **********************************************************

// x^2 and 1/x^2
inline double sq(double x) { return (x * x); }
inline double sqin(double x) { return (1.0 / sq(x)); }

inline double fn_f(double r, double r2m) { return r / (r + r2m); }

inline double fn_beta(double r, double r2m) { return r2m / (r + r2m); }

void param_collect(char **source, int num, map<string, string>& dest) {
  for (int arg = 1; arg < num; ++arg) {
    if (source[arg][0] == '-') {
      dest[source[arg]] = source[arg+1];
    }
  }
}

void param_set(map<string, string>& p_all, map<string, string *>& p_str,
	       map<string, int *>& p_int, map<string, double *>& p_dbl,
	       map<string, bool *>& p_bool) {
  for (pair<string, string> p : p_all) {
    if (p_str.count(p.first)) { *p_str[p.first] = p.second; }
    else if (p_int.count(p.first)) { *p_int[p.first] = atoi(&p.second[0]); }
    else if (p_dbl.count(p.first)) { *p_dbl[p.first] = atof(&p.second[0]); }
    else if (p_bool.count(p.first)) { *p_bool[p.first] = (bool) atoi(&p.second[0]); }
  }
}

// write fields using bbhutil
void write_step(double* fields[], int nfields, char* files[],
		double time, int* shape, int rank, double* coordinates) {
  for (int k = 0; k < nfields; ++k) {
    gft_out_bbox(files[k], time, shape, rank, coordinates, fields[k]);
  }
  return;
}

// get coarsened arrays from fields for writing
void get_write_arr(const vector<double>& field1, const vector<double>& field2,
		    vector<double>& write1, vector<double>& write2,
		    int one_past_last, int savept)
{
  int k, s;
  for (k = 0; k < one_past_last; ++k) {
    s = savept*k;
    write1[k] = field1[s];
    write2[k] = field2[s];
  }
  return;
}

// compute and write mass
void mass_check(const vector<double>& f1, const vector<double>& f2,
		double t, double dr, double rmin, ofstream& out_stream,
		const vector<double>& f_vec, const vector<double>& beta_vec)
{
  double mass = 0.0, rval = rmin;
  int k = 0;
  for (auto val : f2) {
    mass += (0.5*f_vec[k]*(sq(f1[k]) + sq(val)) +
	     beta_vec[k]*f1[k]*val) * 4*M_PI*rval*rval*dr;
    ++k;
    rval += dr;
  }
  out_stream << t <<","<< mass << endl;
  return;
}

// kreiss-oliger dissipation (p.23 choptuik notes)
inline double dissipate(double eps, const vector<double>& u, int ind)
{ return -eps * 0.0625 * ( u[ind+2] -4*u[ind+1] + 6*u[ind]
			   - 4*u[ind-1] + u[ind-2] ); }

// centered differencing operator for fn1*field1 + fn2*field2
// --> mult by 1/(2*dr) to get d(fn1*field1 + fn2*field2)/dr at O(dr^2)
inline double dif_c(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    int ind)
{ return fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1]
    - fn1[ind-1]*field1[ind-1] - fn2[ind-1]*field2[ind-1]; }

// centered differencing operator for r^2*(fn1*field1 + fn2*field2)
// --> mult by 1/(2*dr) to get d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2dif_c(const vector<double>& fn1, const vector<double>& field1,
		      const vector<double>& fn2, const vector<double>& field2,
		      int ind, double dr, double r)
{ return sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
    - sq(r-dr)*(fn1[ind-1]*field1[ind-1] + fn2[ind-1]*field2[ind-1]); }

// forward differencing operator for fn1*field1 + fn2*field2
// --> mult by 1/(2*dr) to get d(fn1*field1 + fn2*field2)/dr at O(dr^2)
inline double dif_f(const vector<double>& fn1, const vector<double>& field1,
		    const vector<double>& fn2, const vector<double>& field2,
		    int ind)
{ return -3*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind])
    + 4*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
    - (fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]); }

// forward differencing operator for r^2*(fn1*field1 + fn2*field2)
// --> mult by 1/(2*dr) to get d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2dif_f(const vector<double>& fn1, const vector<double>& field1,
		      const vector<double>& fn2, const vector<double>& field2,
		      int ind, double dr, double r)
{ return -3*sq(r)*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind])
    + 4*sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
    - sq(r+2*dr)*(fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]); }

// centered differencing operator for r^2*(fn1*field1 + fn2*field2) [using d/d(r^3)]
// --> mult by 1/(2*dr) to get (1/r^2)*d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2d3_c(const vector<double>& fn1, const vector<double>& field1,
		     const vector<double>& fn2, const vector<double>& field2,
		     int ind, double dr, double r)
{ return ( (sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
	  - sq(r-dr)*(fn1[ind-1]*field1[ind-1] + fn2[ind-1]*field2[ind-1]))
	  * 3.0 / (3*r*r + dr*dr) ); }

// forward differencing operator for r^2*(fn1*field1 + fn2*field2) [using d/d(r^3)]
// --> mult by 1/(2*dr) to get (1/r^2)*d(r^2*(fn1*field1 + fn2*field2))/dr at O(dr^2)
inline double r2d3_f(const vector<double>& fn1, const vector<double>& field1,
		     const vector<double>& fn2, const vector<double>& field2,
		     int ind, double dr, double r)
{ return ( (-3*sq(r)*(fn1[ind]*field1[ind] + fn2[ind]*field2[ind])
	    + 4*sq(r+dr)*(fn1[ind+1]*field1[ind+1] + fn2[ind+1]*field2[ind+1])
	    - sq(r+2*dr)*(fn1[ind+2]*field1[ind+2] + fn2[ind+2]*field2[ind+2]))
	   * 3.0 / (3*r*r - 2*dr*dr) ); }

// centered ires_f1 = f1 - ires_c(f1, f2, ind, c1, d1, c2, d2)
//                   - oldf1 - ires_c(oldf1, oldf2, ind, c1, d1, c2, d2)
inline double ires_c(const vector<double>& field1, const vector<double>& field2,
			int ind, double c1, double d1, double c2, double d2)
{ return c1*field1[ind] + c2*field2[ind]
    + d1*(field1[ind+1] - field1[ind-1])
    + d2*(field2[ind+1] - field2[ind-1]); }

// forward ires_f1 = f1 - ires_c(f1, f2, ind, c1, d1, c2, d2)
//                   - oldf1 - ires_c(oldf1, oldf2, ind, c1, d1, c2, d2)
inline double ires_f(const vector<double>& field1, const vector<double>& field2,
			int ind, double c1, double d1, double c2, double d2)
{ return c1*field1[ind] + c2*field2[ind]
    + d1*(-3*field1[ind] + 4*field1[ind+1] - field1[ind+2])
    + d2*(-3*field2[ind] + 4*field2[ind+1] - field2[ind+2]); }

void get_write_arr_ires(const vector<double>& f2, const vector<double>& f3,
			 vector<double>& w2, vector<double>& w3, int last, int savept,
			 const vector<double>& oldf2, const vector<double>& oldf3,
			 vector<double>& ires2, vector<double>& ires3,
			 const vector<double>& f2c1, const vector<double>& f2d1,
			 const vector<double>& f2c2, const vector<double>& f2d2,
			 const vector<double>& f3c1, const vector<double>& f3d1,
			const vector<double>& f3c2, const vector<double>& f3d2)
{
  w2[0] = f2[0];
  w3[0] = f3[0];
  ires2[0] = f2[0] - ires_f(f2, f3, 0, f2c1[0], f2d1[0], f2c2[0], f2d2[0])
    - oldf2[0] - ires_f(oldf2, oldf3, 0, f2c1[0], f2d1[0], f2c2[0], f2d2[0]);
  ires3[0] = f3[0] - ires_f(f3, f2, 0, f3c1[0], f3d1[0], f3c2[0], f3d2[0])
    - oldf3[0] - ires_f(oldf3, oldf2, 0, f3c1[0], f3d1[0], f3c2[0], f3d2[0]);
  int k, s;
  for (k = 1; k < last; ++k) {
    s = savept*k;
    w2[k] = f2[s];
    w3[k] = f3[s];
    ires2[k] = f2[s] - ires_c(f2, f3, s, f2c1[s], f2d1[s], f2c2[s], f2d2[s])
      - oldf2[s] - ires_c(oldf2, oldf3, s, f2c1[s], f2d1[s], f2c2[s], f2d2[s]);
    ires3[k] = f3[s] - ires_c(f3, f2, s, f3c1[s], f3d1[s], f3c2[s], f3d2[s])
      - oldf3[s] - ires_c(oldf3, oldf2, s, f3c1[s], f3d1[s], f3c2[s], f3d2[s]);
  }
  w2[last] = f2[savept*last];
  w3[last] = f3[savept*last];
  return;
}

// **********************************************************
// **********************************************************
//             INITIAL AND BOUNDARY CONDITIONS
// **********************************************************
// **********************************************************

// for gaussian field or sin(coeff*r)*cos(coeff*t)/(coeff*r)
inline double ic_field(double r, double amp, double dsq, double r0)
{ return amp * exp(-(r - r0)*(r - r0)/dsq); }

inline double ic_phi(double r, double amp, double dsq, double r0)
{ return -2 * (r - r0) * amp * exp(-(r - r0)*(r - r0)/dsq) / dsq; }

inline double ic_pi(double r, double amp, double dsq, double r0, double r2m)
{ return ic_phi(r, amp, dsq, r0) + ic_field(r, amp, dsq, r0)/(r*fn_f(r, r2m)); }

inline void phi_bcR(vector<double>& Phi, const vector<double>& oldPhi,
		    int ind, double lambda, double cJ, bool somm_bc) {
  // sommerfeld
  if (somm_bc) {
    Phi[ind] = (lambda*(Phi[ind-1] + oldPhi[ind-1])
		- 0.25*lambda*(Phi[ind-2] + oldPhi[ind-2])
		+ (1 - cJ)*oldPhi[ind]) / (1 + cJ);
  }
  // dirichlet (which means neumann for space derivative)
  else { Phi[ind] = (4*Phi[ind-1] - Phi[ind-2]) / 3.0; } //0; }
  return;
}
inline void pi_bcR(vector<double>& Pi, const vector<double>& oldPi,
		   int ind, double lambda, double cJ, bool somm_bc) {
  // sommerfeld
  if (somm_bc) {
      Pi[ind] = (lambda*(Pi[ind-1] + oldPi[ind-1])
		 - 0.25*lambda*(Pi[ind-2] + oldPi[ind-2])
		 + (1 - cJ)*oldPi[ind]) / (1 + cJ);
  }
  // dirichlet
  else { Pi[ind] = 0; }
  return;
}

// **********************************************************
// **********************************************************
//                         PROGRAM
// **********************************************************
// **********************************************************

int main(int argc, char **argv)
{
  // **********************************************************
  // ******************** PARAMETERS **************************
  // **********************************************************
  
  // user-set parameters
  string outfile = "p1";
  int lastpt = 1000; // grid size
  int save_pt = 1; // write only every (save_pt)th grid point
  int nsteps = 4000; // time steps
  int save_step = 8; // write only every (save_step)th time step
  double lam = 0.25; // dt/dr
  double r2m = 2.0;
  double rmin = 1.5;
  double rmax = 100.0;
  double dspn = 0.5; // dissipation coefficient
  double tol = 0.000000000001; // iterative method tolerance
  int maxit = 25; // max iterations for debugging
  double ic_Dsq = 2.0; // gaussian width
  double ic_r0 = 50.0; // gaussian center
  double ic_Amp = 1.0; // gaussian amplitude
  int check_step = 100; // for monitoring invariant mass
  // note: set bools in command line with integers 1=true or 0=false
  bool zero_pi = false; // zero initial time derivative?
  bool sommerfeld = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  bool write_ires = false; // write ires?
  bool write_res = false; // write res?
  bool write_field = false; // write field?
  bool write_itn = false; // write itn counts?
  // variable to hold constant across resolutions
  string hold_const = "lambda"; // "lambda", "dt", or "dr"

  map<string, string *> p_str {{"-outfile",&outfile},
      {"-hold_const",&hold_const}};
  map<string, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
      {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
      {"-check_step", &check_step}};
  map<string, double *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ic_Dsq",&ic_Dsq},
      {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
  map<string, bool *> p_bool {{"-zero_pi",&zero_pi},
      {"-sommerfeld",&sommerfeld}, {"-dspn_bound",&dspn_bound},
      {"-write_ires",&write_ires}, {"-write_res",&write_res},
      {"-write_field",&write_field}, {"-write_itn",&write_itn}};
  map<string, string> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);

  // check that grid size (lastpt = npts-1) is divisible by save_pt 
  if (lastpt % save_pt != 0) {
    cout << "ERROR: save_pt = " << save_pt << " entered for grid size " << lastpt << endl;
    save_pt -= lastpt % save_pt;
    cout << "--> corrected: save_pt = " << save_pt << endl;
  }
  
  // OBTAIN RESOLUTION FACTORS
  string resn_str;
  int num_resns;
  cout << "enter integer number of resolutions:" << endl;
  cin >> resn_str;
  num_resns = atoi(&resn_str[0]);
  vector<int> resolutions(num_resns);
  for (int k = 0; k < num_resns; ++k) {
    cout << "enter integer resolution factor at level " << k << ":" <<endl;
    cin >> resn_str;
    resolutions[k] = atoi(&resn_str[0]);
  }

  // bbhutil parameters for writing data to sdf
  int lastwrite = lastpt/save_pt;
  int write_shape = lastwrite + 1;
  vector<double> wr_phi(write_shape), wr_pi(write_shape);
  double *field_arr[2] = {&wr_phi[0], &wr_pi[0]};
  int *bbh_shape = &write_shape;
  int bbh_rank = 1;
  double coord_lims[2] = {rmin, rmax};
  double *coords = &coord_lims[0];
  
  int ires_size = ((write_ires) ? write_shape : 1);
  vector<double> iresphi(ires_size, 0.0), irespi(ires_size, 0.0);
  double *ires_arr[2] = {&iresphi[0], &irespi[0]};

  vector<double> field(((write_field) ? write_shape : 1), 0.0);	
  
  // **************************************************************
  // **************************************************************
  //                 LOOP PROGRAM OVER RESOLUTIONS
  // **************************************************************
  // **************************************************************

  int lastpt0 = lastpt; 
  int save_pt0 = save_pt;
  int nsteps0 = nsteps;
  int save_step0 = save_step;
  string outfile0 = outfile;
  double lam0 = lam;
  
  for (int factor : resolutions) {
    if (hold_const == "lambda") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      outfile = to_string(factor) + "-" + outfile0;
    }
    else if (hold_const == "dt") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      lam = lam0 * factor;
      outfile = to_string(factor) + "dr-" + outfile0;
    }
    else if (hold_const == "dr") {
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      lam = lam0 / ((double) factor);
      outfile = to_string(factor) + "dt-" + outfile0;
    }
    else { cout << "ERROR: hold_const must be 'lambda' or 'dt' or 'dr'" << endl; }
  
    // derived parameters
    int npts = lastpt + 1;
    double dr = (rmax - rmin) / ((double) lastpt);
    double dt = lam * dr;
    double somm_coeff = 0.75*lam + 0.5*dt/rmax; // for outer bc
    
  // OUTPUT parameter data
  cout << "\noutfile name = " << outfile << "\ngrid size = " << lastpt << " (" << save_pt
       << "/write)\ntime steps = " << nsteps << " (" << save_step << "/write)\nlambda = "
       << lam << "\nr2m = " << r2m << "\nrmin = " << rmin << "\nrmax = " << rmax
       << "\ndissipation = " << dspn << "\niterative tolerance = " << tol << "\nmaximum iterations = "
       << maxit << "\nic_Dsq = " << ic_Dsq << "\nic_r0 = " << ic_r0 << "\nic_Amp = " << ic_Amp
       << "\nmass check step = " << check_step << "\nmaximum evolution time = " << nsteps*dt
       << "\ndr = " << dr << "\ndt = " << dt << endl;
  ofstream specs;
  string specs_name = outfile + ".txt";
  specs.open(specs_name, ofstream::out);
  specs << "\noutfile name = " << outfile << "\ngrid size = " << lastpt << " (" << save_pt
	<< "/write)\ntime steps = " << nsteps << " (" << save_step << "/write)\nlambda = "
	<< lam << "\nr2m = " << r2m << "\nrmin = " << rmin << "\nrmax = " << rmax
	<< "\ndissipation = " << dspn << "\niterative tolerance = " << tol << "\nmaximum iterations = "
	<< maxit << "\nic_Dsq = " << ic_Dsq << "\nic_r0 = " << ic_r0 << "\nic_Amp = " << ic_Amp
	<< "\nmass check step = " << check_step << "\nmaximum evolution time = " << nsteps*dt
	<< "\ndr = " << dr << "\ndt = " << dt << "\n\noptions:\nzero pi_0 = " << boolalpha << zero_pi
	<< "\nsommerfeld bc = " << sommerfeld << "\ndissipation at bound = " << dspn_bound << endl;
  specs.close();
  
  // **********************************************************
  // ***************** OBJECT DECLARATIONS ********************
  // **********************************************************
  
  // outfiles
  string solname = "sol-" + outfile + ".sdf";
  string outfilePhi_name = "Phi-" + outfile + ".sdf";
  string outfilePi_name = "Pi-" + outfile + ".sdf";
  char *name_arr[2] = {&outfilePhi_name[0], &outfilePi_name[0]};
  string iresPhi_name = "iresPhi-" + outfile + ".sdf";
  string iresPi_name = "iresPi-" + outfile + ".sdf";
  char *iresname_arr[2] = {&iresPhi_name[0], &iresPi_name[0]};  
  string mass_file = "mass-" + outfile + ".csv";
  ofstream ofs_mass;
  ofs_mass.open(mass_file, ofstream::out);
  ofs_mass << "save=" << save_step <<","<< "check=" << check_step << endl;
  string itn_file = "itns-" + outfile + ".csv";
  ofstream ofs_itn;
  if (write_itn) { ofs_itn.open(itn_file, ofstream::out); }

  // fields and residuals
  vector<double> phi(npts, 0.0), pi(npts, 0.0);
  vector<double> old_phi(npts, 0.0), old_pi(npts, 0.0);
  vector<double> bphi(npts, 0.0), bpi(npts, 0.0);
  vector<double> resphi(npts, 0.0), respi(npts, 0.0);
  vector<double> f(npts, 0.0), beta(npts, 0.0);
  ires_size = ((write_ires) ? npts : 1);
  vector<double> phic1(ires_size, 0.0),  phic2(ires_size, 0.0),
    pic1(ires_size, 0.0), pic2(ires_size, 0.0),
    d1(ires_size, 0.0), d2(ires_size, 0.0);

  // *********************************************
  // **************** DEBUG **********************
  // *********************************************
  string resphi_fname = "resPhi-" + outfile + ".sdf";
  string respi_fname = "resPi-" + outfile + ".sdf";
  char *resname_arr[2] = {&resphi_fname[0], &respi_fname[0]};
  int maxit_count = 0;

  time_t start_time = time(NULL); // time for rough performance measure
  
  // **********************************************************
  // ******************* INITIAL DATA ************************
  // **********************************************************

  int i, j, itn = 0; // declare loop integers
  double res = 1.0; // declare residual indicator
  double r = rmin, t = 0.0; // declare position and time variables
  for (j = 0; j < npts; ++j) {
    phi[j] = ic_phi(r, ic_Amp, ic_Dsq, ic_r0);
    if (!zero_pi) { pi[j] = ic_pi(r, ic_Amp, ic_Dsq, ic_r0, r2m); }
    f[j] = fn_f(r, r2m);
    beta[j] = fn_beta(r, r2m);
    if (write_ires) {
      phic1[j] = -0.5*dt*sq(beta[j])/r2m;
      phic2[j] = -phic1[j];
      pic1[j] = phic1[j] + dt*beta[j]/r;
      pic2[j] = phic2[j] + dt*f[j]/r;
      d1[j] = 0.25*lam*beta[j];
      d2[j] = 0.25*lam*f[j];
    }
    if ((write_field) && (j%save_pt == 0)) {
      field[j/save_pt] = ic_field(r, ic_Amp, ic_Dsq, ic_r0); }
    r += dr;
  }

  // **********************************************************
  // ******************* TIME STEPPING ************************
  // *******************   & WRITING   ************************
  // **********************************************************
  
  gft_set_multi(); // start bbhutil file i/o
  for (i = 0; i < nsteps; ++i) {
    t = i * dt;
    // *************  WRITE the fields (& res, ires) at t(n)  **************
    if (i % save_step == 0) {
      if (write_ires) {
	get_write_arr_ires(phi, pi, wr_phi, wr_pi, lastwrite, save_pt,
			   old_phi, old_pi, iresphi, irespi, phic1, d1,
			   phic2, d2, pic1, d1, pic2, d2);
	write_step(ires_arr, 2, iresname_arr, t, bbh_shape, bbh_rank, coords);
      }
      else { get_write_arr(phi, pi, wr_phi, wr_pi, lastwrite+1, save_pt); }
      
      write_step(field_arr, 2, name_arr, t, bbh_shape, bbh_rank, coords);
      
      if (write_res) {
	get_write_arr(resphi, respi, wr_phi, wr_pi, lastwrite+1, save_pt);
	write_step(field_arr, 2, resname_arr, t, bbh_shape, bbh_rank, coords);
      }
      if (write_field) { gft_out_bbox(&solname[0], t, bbh_shape, bbh_rank,
				      coords, &field[0]); }	
    }
    // now set old_phi/pi to t(n) so that phi/pi can be updated to t(n+1)
    old_phi = phi;
    old_pi = pi;

    // **********************************************************
    //         SOLVE EVOLUTION EQUATION Ax = b for x(n+1)
    // **********************************************************   
    // create rhs vec b in A_lhs.x(n+1) = A_rhs.x(n) := b
    r = rmin;
    bphi[0] = old_phi[0] + 0.25*lam*dif_f(beta, old_phi, f, old_pi, 0);
    bpi[0] = old_pi[0] + 0.25*lam*sqin(r)*r2dif_f(beta, old_pi, f, old_phi, 0, dr, r);
    for (j = 1; j < lastpt; ++j) {
      r += dr;
      bphi[j] = old_phi[j] + 0.25*lam*dif_c(beta, old_phi, f, old_pi, j);
      bpi[j] = old_pi[j] + 0.25*lam*sqin(r)*r2dif_c(beta, old_pi, f, old_phi, j, dr, r);
    }
    // reset itn and set res > tol to enter GAUSS-SEIDEL ITERATIVE SOLVER
    itn = 0, res = 1.0;
    while (res > tol) {
      
      r = rmin;
      // UPDATE INTERIOR and collect residuals
      phi[0] = bphi[0] + 0.25*lam*dif_f(beta, phi, f, pi, 0);
      pi[0] = bpi[0] + 0.25*lam*sqin(r)*r2dif_f(beta, pi, f, phi, 0, dr, r);
      r += dr;
      phi[1] = bphi[1] + 0.25*lam*dif_c(beta, phi, f, pi, 1);
      pi[1] = bpi[1] + 0.25*lam*sqin(r)*r2dif_c(beta, pi, f, phi, 1, dr, r);
      for (j = 2; j < lastpt; ++j) {
        r += dr;
	phi[j] = bphi[j] + 0.25*lam*dif_c(beta, phi, f, pi, j);
	pi[j] = bpi[j] + 0.25*lam*sqin(r)*r2dif_c(beta, pi, f, phi, j, dr, r);
	resphi[j-1] = abs(phi[j-1] - bphi[j-1] -
			  0.25*lam*dif_c(beta, phi, f, pi, j-1));
	respi[j-1] = abs(pi[j-1] - bpi[j-1] - sqin(r-dr)*
			 0.25*lam*r2dif_c(beta, pi, f, phi, j-1, dr, r-dr));
      }
      resphi[0] = abs(phi[0] - bphi[0] -
		      0.25*lam*dif_f(beta, phi, f, pi, 0));
      respi[0] = abs(pi[0] - bpi[0] - sqin(rmin)*
		     0.25*lam*r2dif_f(beta, pi, f, phi, 0, dr, rmin));
      // UPDATE BOUNDARY
      phi_bcR(phi, old_phi, lastpt, lam, somm_coeff, sommerfeld);
      pi_bcR(pi, old_pi, lastpt, lam, somm_coeff, sommerfeld);
      
      resphi[lastpt-1] = abs(phi[lastpt-1] - bphi[lastpt-1] -
			     0.25*lam*dif_c(beta, phi, f, pi, lastpt-1));
      respi[lastpt-1] = abs(pi[lastpt-1] - bpi[lastpt-1] - sqin(r)*
			    0.25*lam*r2dif_c(beta, pi, f, phi, lastpt-1, dr, r));
      
      // CHECK RESIDUAL
      res = max(*max_element(resphi.begin(), resphi.end()), *max_element(respi.begin(), respi.end())); // can also use 1-norm or 2-norm

      ++itn; 
      if (itn % maxit == 0) {
	res = 0.0;
	++maxit_count;
	if (i % 500*factor == 0) { cout << i << " res= " << res << " at " << itn << endl; }
      }
    }   
    if (write_itn) { ofs_itn << itn << endl; } // record itn count
    // ****************** ITERATIVE SOLUTION COMPLETE ******************
    
    // ****************** kreiss-oliger DISSIPATION ********************
    // at ind next to boundaries can call dissipate on ind+/-1 or ind+/-2
    if (dspn_bound) {
      phi[0] += dissipate(dspn, old_phi, 2);
      pi[0] += dissipate(dspn, old_pi, 2);
      phi[1] += dissipate(dspn, old_phi, 2);
      pi[1] += dissipate(dspn, old_pi, 2);
    }
    for (j = 2; j < lastpt-1; ++j) {
      phi[j] += dissipate(dspn, old_phi, j);
      pi[j] += dissipate(dspn, old_pi, j);
    }

    // ****************** WRITE MASS & update field **********************
    if (i % check_step*save_step == 0) { mass_check(phi, pi, t, dr, rmin, ofs_mass, f, beta); }
    if (write_field) {
      int k, s;
      for (k = 0; k < write_shape; ++k) {
	s = k*save_pt;
	field[k] += dt * (beta[s]*phi[s] + f[s]*pi[s]);
      }
    } 
  }
  // ******************** DONE TIME STEPPING *********************
  
  // write final time step
  if (nsteps % save_step == 0) {
    get_write_arr(phi, pi, wr_phi, wr_pi, lastwrite+1, save_pt);
    write_step(field_arr, 2, name_arr, nsteps*dt, bbh_shape, bbh_rank, coords);
  }
  // close outfiles
  gft_close_all();
  ofs_mass.close();
  if (write_itn) { ofs_itn.close(); }
  // print resolution runtime
  cout << difftime(time(NULL), start_time) << " seconds elapsed" << endl;
  //*************DEBUG************
  cout << maxit_count << " steps reached maxit=" << maxit << "\n" << endl;
  
  }
  // ******************** DONE LOOPING OVER RESOLUTIONS *********************
  
  return 0;
}
