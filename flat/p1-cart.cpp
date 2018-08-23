/* 

radial wave equation for massless scalar field in flat space

compile with:

   g++ -std=c++11 -g -Wall -O2 p1-flat.cpp -o p1-flat -lbbhutil

input parameters in terminal as:
(do NOT include .sdf extension in outfile_name)

./p1-flat <lastpt> <save_pt> <nsteps> <save_step> <lam> <rmax> <outfile_name>
          <iterative_tolerance> <dissipation_coeff> <steps_per_mass_check>
	  <ic_Dsq> <ic_Amp> <ic_r0>
*/

#include <iostream>
#include <algorithm> // for max_element()
#include <vector> // for everything
#include <fstream> // for mass/itn files
#include <ctime> // for quick time analysis
#include <string> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include <cmath> // for ICs
#include "bbhutil.h" // for output to .sdf

#define _USE_MATH_DEFINES
#ifndef M_PI
cerr << "NO PI\n";
const double M_PI = 4.0*atan(1.0);
#endif


using namespace std;

// **********************************************************
//                       FUNCTIONS
// **********************************************************

// BOUNDARY STENCIL FOR NO BC
  /*
inline void phi_bound0sten(vector<double>& Phi, const vector<double>& Pi, const vector<double>& oldPhi, const vector<double>& oldPi, double lambda) {
  Phi[0] = oldPhi[0] - 0.25*lambda*(3*oldPi[0] - 4*oldPi[1] + oldPi[2] + 3*Pi[0] - 4*Pi[1] + Pi[2]);
  return;
}
inline void pi_bound0sten(const vector<double>& Phi, vector<double>& Pi, const vector<double>& oldPhi, const vector<double>& oldPi, double lambda){
  Pi[0] = oldPi[0] - 0.25*lambda*(3*oldPhi[0] - 4*oldPhi[1] + oldPhi[2] + 3*Phi[0] - 4*Phi[1] + Phi[2]);
  return;
}
inline void phi_boundRsten(vector<double>& Phi, const vector<double>& Pi, const vector<double>& oldPhi, const vector<double>& oldPi, double lambda, int ind){
  Phi[ind] = oldPhi[ind] + 0.25*lambda*(3*oldPi[ind] - 4*oldPi[ind-1] + oldPi[ind-2] + 3*Pi[ind] - 4*Pi[ind-1] + Pi[ind-2]);
  return;
}
inline void pi_boundRsten(const vector<double>& Phi, vector<double>& Pi, const vector<double>& oldPhi, const vector<double>& oldPi, double lambda, int ind){
  Pi[ind] = oldPi[ind] + 0.25*lambda*(3*oldPhi[ind] - 4*oldPhi[ind-1] + oldPhi[ind-2] + 3*Phi[ind] - 4*Phi[ind-1] + Phi[ind-2]);
  return;
}
  */

// write fields using bbhutil
void write_step(double* fields[], int nfields, char* files[], double time, int* shape, int rank, double* coordinates) {
  for (int k = 0; k < nfields; ++k) {
    gft_out_bbox(files[k], time, shape, rank, coordinates, fields[k]);
  }
  return;
}

// get coarsened arrays from field for writing
void get_write_arrs(const vector<double>& field1, const vector<double>& field2, const vector<double>& field3,
		    vector<double>& write1, vector<double>& write2, vector<double>& write3, int one_past_last, int savept)
{
  for (int k = 0; k < one_past_last; ++k) {
    write1[k] = field1[savept*k];
    write2[k] = field2[savept*k];
    write3[k] = field3[savept*k];
  }
  return;
}

// get coarsened arrays from field while using previous time step to start ires computation
void get_write_arrs_ires(const vector<double>& field1, const vector<double>& field2, const vector<double>& field3,
			 vector<double>& write1, vector<double>& write2, vector<double>& write3, int last, int savept,
			 const vector<double>& oldfield2, const vector<double>& oldfield3,
			 vector<double>& ires2, vector<double>& ires3, double lam)
{
  write1[0] = field1[0];
  write2[0] = field2[0];
  write3[0] = field3[0];
  for (int k = 1; k < last; ++k) {
    write1[k] = field1[savept*k];
    write2[k] = field2[savept*k];
    write3[k] = field3[savept*k];
    // now old_phi/pi is from t_(n-1), so can get indep. residual
    // from scheme with standard 2nd-order centered scheme in both
    // space and time: D_t.x = D_r.y   (+1/r term for Pi) -->
    // x_(n+1)[j] - x_(n-1)[j] + lam*(y_(n)[j+1] - y_(n)[j-1] = 0
    ires2[k] = oldfield2[savept*k] - lam*field3[savept*(k-1)] + lam*field3[savept*(k+1)];
    ires3[k] = oldfield3[savept*k] - lam*field2[savept*(k-1)] + 4*lam*field2[k]/((double) k) + lam*field2[savept*(k+1)];
    // we will then set iresphi/pi = phi/pi - iresphi/pi
    // to obtain the residual of the independent FDA equation above
    // once we have updated phi/pi to t(n+1)
  }
  write1[last] = field1[savept*last];
  write2[last] = field2[savept*last];
  write3[last] = field3[savept*last];
  return;
}

void finish_ires_comp(vector<double>& ires1, vector<double>& ires2,
		      const vector<double>& field1, const vector<double>& field2, int last, int savept) {
  for (int k = 1; k < last; ++k) {
    ires1[k] = field1[savept*k] - ires1[k];
    ires2[k] = field2[savept*k] - ires2[k];
  }
  return;
}

// update field using auxiliary fields
void update_field(vector<double>& field1, const vector<double>& field2, const vector<double>& field3, double dt) {
  // for flat just fied[j] += pi[j]*dt
  int k = 0;
  for (auto val : field3) { field1[k++] += val * dt; }
}

// update field using auxiliary fields while writing mass
void update_field_mass_check(vector<double>& field1, const vector<double>& field2, const vector<double>& field3,
			     double dt, double dr, ofstream& out_stream)
{
  // for flat just fied[j] += pi[j]*dt
  double mass = 0.0;
  int k = 0;
  for (auto val : field3) {
    mass += 2 * M_PI * (k*dr)*(k*dr) * (field2[k]*field2[k] + field3[k]*field3[k]) * dr;
    field1[k++] += val * dt;
  }
  out_stream << mass << endl;
}

// 2-norm of vector
double norm2(const vector<double>& vec) {
  double norm = 0.0;
  for (auto x : vec) { norm += x * x; }
  return sqrt(norm);
}

// dissipation term from p.23 choptuik notes
inline double dissipate(double eps, const vector<double>& u, int ind)
{ return -eps * 0.0625 * ( u[ind+2] -4*u[ind+1] + 6*u[ind]
			   - 4*u[ind-1] + u[ind-2] ); }

// **********************************************************
// **********************************************************
//             INITIAL AND BOUNDARY CONDITIONS
// **********************************************************
// **********************************************************

// for field sin(coeff*r)*cos(coeff*t)/(coeff*r)
inline double ic_field(double r, double amp, double dsq, double r0) {
  if (r == 0.0) { return 0.0; }
  //else { return sin(coeff*r) / (coeff*r); }
  
  else { return amp * exp(-(r - r0)*(r - r0) / dsq); }
  //double coeff = M_PI / r0;
  //return sin(coeff*r);
}
inline double ic_phi(double r, double amp, double dsq, double r0) {
  if (r == 0.0) { return 0.0; }
  //else { return (cos(coeff*r) - (sin(coeff*r)/(coeff*r))) / r; }
  else { return -2 * (r - r0) * amp * exp(-(r - r0)*(r - r0) / dsq) / dsq; }
  //double coeff = M_PI / r0;
  //return coeff * cos(coeff*r);
}
inline double ic_pi(double r, double amp, double dsq, double r0) {
  if (r == 0.0) { return 0.0; }
  else { return amp * exp(-(r - r0)*(r - r0) / dsq) * (-2*(r - r0)/dsq /*+ (1/r)*/); }
  //return 0.0;
}

inline void phi_bc0(vector<double>& Phi, double lambda) {
  Phi[0] = 0.0;//(4*Phi[1] - Phi[2]) / 3.0;
  return;
}

inline void phi_bcR(vector<double>& Phi, const vector<double>& oldPhi, int ind, double lambda) {
  //Phi[ind] = 0.0;//(4*Phi[ind-1] - Phi[ind-2]) / 3.0;//cos(coeff*t)/rmax;
  double cJ = lambda * ( 0.75 /*+ 0.5/((double) ind)*/ );
  Phi[ind] = (lambda*(Phi[ind-1] + oldPhi[ind-1]) - 0.25*lambda*(Phi[ind-2] + oldPhi[ind-2])
	      + (1 - cJ)*oldPhi[ind]) / (1 + cJ);
  return;
}
inline void pi_bc0(vector<double>& Pi, double lambda) {
  Pi[0] = (4*Pi[1] - Pi[2]) / 3.0;//-coeff*sin(coeff*t);
  return;
}
inline void pi_bcR(vector<double>& Pi, const vector<double>& oldPi, int ind, double lambda) {
  double cJ = lambda * ( 0.75 /*+ 0.5/((double) ind)*/ );
  Pi[ind] = (lambda*(Pi[ind-1] + oldPi[ind-1]) - 0.25*lambda*(Pi[ind-2] + oldPi[ind-2])
	     + (1 - cJ)*oldPi[ind]) / (1 + cJ);
  return;
}

inline double cf(double lambda) { return 0.25 * lambda; }

inline double cG(double r, double delta_t) { return 0.0; }//delta_t / r; }

// **********************************************************
//                       PROGRAM
// **********************************************************

int main(int argc, char **argv)
{
  int maxit = 20; // max iterations for debugging

  // **********************************************************
  //                       PARAMETERS
  // **********************************************************
  // user-set parameters
  int lastpt = 2000; // grid size
  int save_pt = 2; // write only every (save_pt)th grid point
  int nsteps = 4000; // time steps
  int save_step = 4; // write only every (save_step)th time step
  double lam = 0.5; // dt/dr
  double rmax = 100.0;
  string outfile_name = "flat";
  double tol = 0.000000000001; // iterative method tolerance
  double dspn = 0.0; // dissipation coefficient
  int check_step = 100; // for monitoring invariant mass
  double ic_Dsq = 2.0; // gaussian width
  double ic_Amp = 1.0; // gaussian amplitude
  double ic_r0 = 50.0; // gaussian center
  if (argc > 7) {
    lastpt = atoi(&argv[1][0]);
    save_pt = atoi(&argv[2][0]);
    nsteps = atoi(&argv[3][0]);
    save_step = atoi(&argv[4][0]);
    lam = atof(&argv[5][0]);
    rmax = atof(&argv[6][0]);
    outfile_name = argv[7];
    if (argc > 10) {
      tol = atof(&argv[8][0]);
      dspn = atof(&argv[9][0]);
      check_step = atoi(&argv[10][0]);
      if (argc == 14) {
	ic_Dsq = atof(&argv[11][0]);
	ic_Amp = atof(&argv[12][0]);
	ic_r0 = atof(&argv[13][0]);
	}
    }
  }
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
  string hold_const;
  cout << "enter variable to hold constant (lambda, dt, dr):" << endl;
  cin >> hold_const;

  // ask about saving independent residuals and itn counts
  bool write_ires = false;
  string ires_str;
  cout << "want to save indpendent residuals? (1/0)" << endl;
  cin >> ires_str;
  write_ires = (bool) atoi(&ires_str[0]);
  bool write_itn = false;
  string itn_str;
  cout << "want to save iteration counts? (1/0)" << endl;
  cin >> itn_str;
  write_itn = (bool) atoi(&itn_str[0]);

  int lastpt0 = lastpt; 
  int save_pt0 = save_pt;
  int nsteps0 = nsteps;
  int save_step0 = save_step;
  string outfile_name0 = outfile_name;
  double lam0 = lam;
  
    //***************************************************************
    //              LOOP PROGRAM OVER RESOLUTIONS
    // **************************************************************
  
  for (int factor : resolutions) {
    if (hold_const == "lambda") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      outfile_name = outfile_name0 + to_string(factor);
    }
    else if (hold_const == "dt") {
      lastpt = lastpt0 * factor;
      save_pt = save_pt0 * factor;
      lam = lam0 * factor;
      outfile_name = outfile_name0 + "-dr" + to_string(factor);
    }
    else if (hold_const == "dr") {
      nsteps = nsteps0 * factor;
      save_step = save_step0 * factor;
      lam = lam0 / ((double) factor);
      outfile_name = outfile_name0 + "-dt" + to_string(factor);
    }
    else { cout << "ERROR: hold_const must be 'lambda' or 'dt' or 'dr'" << endl; }
  

    int lastwrite = lastpt/save_pt;
    int write_shape = lastwrite + 1;
    int npts = lastpt + 1;
    double dr = rmax / ((double) lastpt);
    double dt = lam * dr;
  
  // OUTPUT parameter data
  cout << "\ngrid size = " << lastpt << "   ->   dr = " << dr << endl;
  cout << "grid points per write = " << save_pt << endl;
  cout << "time steps = " << nsteps << "   ->   dt = " << dt << endl;
  cout << "time steps per write = " << save_step << endl;
  cout << "lambda = " << lam << endl;
  cout << "max r = " << rmax << ",\tmax t = " << nsteps*dt << endl;
  cout << "outfile name = " << outfile_name << endl;
  cout << "iterative tolerance = " << tol << endl;
  cout << "dissipation coefficient = " << dspn << endl;
  cout << "time steps per mass check = " << check_step << "\n" << endl;
  
  // **********************************************************
  //                   OBJECT DECLARATIONS
  // **********************************************************
  
  // output file names
  string outfilePhi_name = outfile_name;
  string outfilePi_name = outfile_name;
  string iresPhi_name = outfile_name;
  string iresPi_name = outfile_name;
  string mass_file = outfile_name;
  string itn_file = outfile_name;
  outfile_name += "-sol.sdf";
  outfilePhi_name += "-Phi.sdf";
  outfilePi_name += "-Pi.sdf";
  iresPhi_name += "-PhiRes.sdf";
  iresPi_name += "-PiRes.sdf";
  mass_file += "-mass.csv";
  itn_file += "-itns.csv";

  ofstream ofs_mass;
  ofs_mass.open(mass_file, ofstream::out);
  ofs_mass << check_step << "\n" << endl;
  ofstream ofs_itn;
  if (write_itn) { ofs_itn.open(itn_file, ofstream::out); }

  // fields and residuals
  vector<double> field(npts, 0.0);
  vector<double> phi(npts, 0.0), pi(npts, 0.0);
  vector<double> old_phi(npts, 0.0), old_pi(npts, 0.0);
  vector<double> bphi(npts, 0.0), bpi(npts, 0.0);
  vector<double> resphi(npts, 0.0), respi(npts, 0.0);
  double res = 1.0;
  int itn = 0;
  vector<double> iresphi(write_shape, 0.0), irespi(write_shape, 0.0);
  vector<double> write_field(write_shape);
  vector<double> write_phi(write_shape), write_pi(write_shape);

  // bbhutil parameters for writing data to sdf
  char *outfile = &outfile_name[0];
  char *outfilePhi = &outfilePhi_name[0];
  char *outfilePi = &outfilePi_name[0];
  char *iresphi_out = &iresPhi_name[0];
  char *irespi_out = &iresPi_name[0];
  double *field_arr[3] = {&write_field[0], &write_phi[0], &write_pi[0]};
  char *name_arr[3] = {outfile, outfilePhi, outfilePi};
  double *ires_arr[2] = {&iresphi[0], &irespi[0]};
  char *iresname_arr[2] = {iresphi_out, irespi_out};
  int bbh_rank = 1;
  int *bbh_shape = &write_shape;
  double coord_lims[2] = {0.0, rmax};
  double *coords = &coord_lims[0];
  gft_set_multi();

  // *********************************************
  // **************** DEBUG **********************
  // *********************************************
  string resphi_fname = "A-" + iresPhi_name;
  string respi_fname = "A-" + iresPi_name;
  char *resphi_out = &resphi_fname[0];
  char *respi_out = &respi_fname[0];
  double *res_arr[2] = {&resphi[0], &respi[0]};
  char *resname_arr[2] = {resphi_out, respi_out};

  time_t start_time = time(NULL); // time for rough performance measure

  // **********************************************************
  // **********************************************************
  //                       TIME STEPPING
  // **********************************************************
  // **********************************************************
  
  // *** INITIAL DATA ***
  int i, j;
  // wavenumber: double ic_k = 2 * M_PI / rmax;
  for (i = 0; i < npts; ++i) {
    field[i] = ic_field(i*dr, ic_Amp, ic_Dsq, ic_r0);
    phi[i] = ic_phi(i*dr, ic_Amp, ic_Dsq, ic_r0);
    pi[i] = ic_pi(i*dr, ic_Amp, ic_Dsq, ic_r0);
  }

  // *** EVOLUTION (w=1 s.o.r.) ***
  double t;
  for (i = 0; i < nsteps; ++i) {
    t = i * dt;
    // *************  WRITE the fields at t(n)  **************
    // *************  & start ires computation  **************
    if (i % save_step == 0) {
      if (write_ires && i > 0) {
	get_write_arrs_ires(field, phi, pi, write_field, write_phi, write_pi,
			    lastwrite, save_pt, old_phi, old_pi, iresphi, irespi, lam);
      }
      else { get_write_arrs(field, phi, pi, write_field, write_phi, write_pi, lastwrite+1, save_pt); }
      write_step(field_arr, 3, name_arr, t, bbh_shape, bbh_rank, coords);
    }

    // now set old_phi/pi to t(n) so that phi/pi can be updated to t(n+1)
    old_phi = phi;
    old_pi = pi;

    // **********************************************************
    // **********************************************************
    //         SOLVE EVOLUTION EQUATION Ax = b for x(n+1)
    // **********************************************************
    // **********************************************************
    
    // create rhs vec b in A_lhs.x(n+1) = A_rhs.x(n) := b (don't need b at boundaries)
    for (j = 1; j < lastpt; ++j) {
      bphi[j] = old_phi[j] - cf(lam)*old_pi[j-1] + cf(lam)*old_pi[j+1];
      bpi[j] = old_pi[j] - cf(lam)*old_phi[j-1] + cG(j*dr, dt)*old_phi[j] + cf(lam)*old_phi[j+1];
    }

    // reset itn and set res > tol to enter GAUSS-SEIDEL ITERATIVE SOLVER (Ax = b +/- tol)
    itn = 0, res = 1.0;
    while (res > tol && itn < maxit) {
      
      // UPDATE BOUNDARY
      phi_bc0(phi, lam);
      pi_bc0(pi, lam);
      double phi_bound0 = phi[0];
      double pi_bound0 = pi[0];
      
      // UPDATE INTERIOR and collect residuals
      phi[1] = bphi[1] - cf(lam)*pi[0] + cf(lam)*pi[2];
      pi[1] = bpi[1] - cf(lam)*phi[0] + cG(dr, dt)*phi[1] + cf(lam)*phi[2];
      for (j = 2; j < lastpt; ++j) {
	phi[j] = bphi[j] - cf(lam)*pi[j-1] + cf(lam)*pi[j+1];
	pi[j] = bpi[j] - cf(lam)*phi[j-1] + cG(j*dr, dt)*phi[j] + cf(lam)*phi[j+1];
	resphi[j-1] = abs(bphi[j-1] - phi[j-1] - cf(lam)*pi[j-2] + cf(lam)*pi[j]);
	respi[j-1] = abs(bpi[j-1] - pi[j-1] - cf(lam)*phi[j-2] + cG((j-1)*dr, dt)*phi[j-1] + cf(lam)*phi[j]);
      }
      
      // UPDATE BOUNDARY and collect final residuals
      phi_bc0(phi, lam);
      pi_bc0(pi, lam);
      phi_bcR(phi, old_phi, lastpt, lam);
      pi_bcR(pi, old_pi, lastpt, lam);

      resphi[0] = abs(phi[0] - phi_bound0);
      respi[0] = abs(pi[0] - pi_bound0);
      resphi[lastpt-1] = abs(bphi[lastpt-1] - phi[lastpt-1] - cf(lam)*pi[lastpt-2] + cf(lam)*pi[lastpt]);
      respi[lastpt-1] = abs(bpi[lastpt-1] - pi[lastpt-1] - cf(lam)*phi[lastpt-2] + cG(rmax-dr, dt)*phi[lastpt-1] + cf(lam)*phi[lastpt]);
      
      // CHECK RESIDUAL
      res = max(*max_element(resphi.begin(), resphi.end()), *max_element(respi.begin(), respi.end()));
      //res = max(norm2(resphi), norm2(respi));

      ++itn;
      
      // TESTING ********* PRINT RESIDUAL IF ITERATION PROLONGED
      if (itn % maxit == 0) { cout << i << " res= " << res << " at " << itn << endl; }
      
      // TESTING ************ DEBUG ********** WRITE RESIDUAL
      //write_step(res_arr, 2, resname_arr, t + itn*dt, bbh_shape, bbh_rank, coords);
    }
    
    if (write_itn) { ofs_itn << itn << endl; } // record itn count

    // ****************** ITERATIVE SOLUTION COMPLETE ******************
    
    // ****************** kreiss-oliger DISSIPATION ********************
    // at ind next to boundaries can call dissipate on ind+/-1 or ind+/-2
    phi[0] += dissipate(dspn, old_phi, 2);
    pi[0] += dissipate(dspn, old_pi, 2);
    phi[1] += dissipate(dspn, old_phi, 2);
    pi[1] += dissipate(dspn, old_pi, 2);
    for (j = 2; j < lastpt-1; ++j) {
      phi[j] += dissipate(dspn, old_phi, j);
      pi[j] += dissipate(dspn, old_pi, j);
    }
    phi[lastpt-1] += dissipate(dspn, old_phi, lastpt-2);
    pi[lastpt-1] += dissipate(dspn, old_pi, lastpt-2);
    phi[lastpt] += dissipate(dspn, old_phi, lastpt-2);
    pi[lastpt] += dissipate(dspn, old_pi, lastpt-2);
    
    // APPLY BOUNDARY CONDITIONS
    phi_bc0(phi, lam);
    pi_bc0(pi, lam);
    phi_bcR(phi, old_phi, lastpt, lam);
    pi_bcR(pi, old_pi, lastpt, lam);

    // ****************** SOLUTION DISSIPATION COMPLETE ******************

    // ****************** MORE WRITING **********************
    // also update field by integrating next step of auxiliary fields
    if (i % save_step == 0) {
      if (write_ires && i > 0) {
	// write ires (currently not doing ires at boundaries)
        finish_ires_comp(iresphi, irespi, phi, pi, lastwrite, save_pt);
	write_step(ires_arr, 2, iresname_arr, t + dt, bbh_shape, bbh_rank, coords);
      }

      // write mass each save_step*check_step steps
      if (i % check_step*save_step == 0) { update_field_mass_check(field, phi, pi, dt, dr, ofs_mass);}
      else { update_field(field, phi, pi, dt); }
    }  
    else { update_field(field, phi, pi, dt); }
  }
  
  // write final time step
  if (nsteps % save_step == 0) {
    get_write_arrs(field, phi, pi, write_field, write_phi, write_pi, lastwrite+1, save_pt);
    write_step(field_arr, 3, name_arr, nsteps*dt, bbh_shape, bbh_rank, coords);
  }
  // *************** WRITING COMPLETE **************
  // close output files
  gft_close_all();
  ofs_mass.close();
  if (write_itn) { ofs_itn.close(); }

  cout << difftime(time(NULL), start_time) << " seconds elapsed" << endl;
  
  }
  
  return 0;
}
