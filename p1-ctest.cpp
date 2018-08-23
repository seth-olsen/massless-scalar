/*

this program will compute

q(t) = ||u_4h(x,t) - u_2h(x,t)||/||u_2h(x,t) - u_h(x,t)||

with (u_4h, u_2h, u_h) from resolution levels (0, 1, 2) 
and save to outfile (.csv) with the line format:
   step_number (= t/dt) , q_phi(t) , q_pi(t) , q_phi+pi(t)

input parameters in terminal as:
(do NOT include .sdf extension in outfile)

 ./p1-ctest <outfile> <lastpt> <save_pt> <nsteps> <save_step>
 <lam> <r2m> <rmin> <rmax> <dspn> <tol> <maxit> <ic_Dsq> 
 <ic_r0> <ic_Amp> <check_step> <zero_pi> <sommerfeld> 
 <dspn_bound> <write_ires> <write_itn> <hold_const>
 <same_times> <same_grids>

where the value of any <parameter> can be set with the
following ordered pair of command line arguments:

 -parameter parameter_value

default values can be found at the start of main()
*/

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cmath>
#include <fstream>
#include "bbhutil.h"
using namespace std;

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

// read fields using bbhutil
void read_step(char *files[], int times[], double *fields[], int nfields) {
  for (int k = 0; k < nfields; ++k) {
    gft_read_brief(files[k], times[k], fields[k]);
  }
  return;
}

int main(int argc, char **argv)
{
  // coarse simulation parameters
  string outfile = "p1";
  string pre1 = "Phi-", pre2 = "Pi-";
  int lastpt = 1000; // grid size
  int save_pt = 1; // write only every (save_pt)th grid point
  int nsteps = 4000; // time steps
  int save_step = 8; // write only every (save_step)th time step
  double lam = 0.25; // dt/dr
  double r2m = 2.0;
  double rmin = 1.5;
  double rmax = 100.0;
  double dspn = 0.8; // dissipation coefficient
  double tol = 0.000000000001; // iterative method tolerance
  int maxit = 25; // max iterations for debugging
  double ic_Dsq = 2.0; // gaussian width
  double ic_r0 = 50.0; // gaussian center
  double ic_Amp = 1.0; // gaussian amplitude
  bool zero_pi = false; // zero initial time derivative?
  bool sommerfeld = true; // sommerfeld condition at outer bound?
  bool dspn_bound = false; // dissipate boundary points?
  // variable to hold constant across resolutions
  string hold_const = "lambda"; // "lambda", "dt", or "dr"
  bool same_times = true;
  bool same_grids = true;
  // resolution factors
  int resn0 = 8, resn1 = 16, resn2 = 32; // in order of priority
  int *resns[3] = {&resn0, &resn1, &resn2};

  // get parameters from command line
  map<string, string *> p_str {{"-outfile",&outfile}, {"-pre1",&pre1},
      {"-pre2",&pre2}, {"-hold_const",&hold_const}};
  map<string, int *> p_int {{"-lastpt",&lastpt}, {"-save_pt", &save_pt},
      {"-nsteps", &nsteps}, {"-save_step",&save_step}, {"-maxit",&maxit},
      {"-resn0", &resn0}, {"-resn1", &resn1}, {"-resn2", &resn2}};
  map<string, double *> p_dbl {{"-lam",&lam}, {"-r2m",&r2m}, {"-rmin",&rmin},
      {"-rmax",&rmax}, {"-dspn",&dspn}, {"-tol",&tol}, {"-ic_Dsq",&ic_Dsq},
      {"-ic_r0",&ic_r0}, {"-ic_Amp",&ic_Amp}};
  map<string, bool *> p_bool {
    {"-sommerfeld",&sommerfeld}, {"-dspn_bound",&dspn_bound},
    {"-same_times",&same_times}, {"-same_grids",&same_grids}};
  map<string, string> params;
  param_collect(argv, argc, params);
  param_set(params, p_str, p_int, p_dbl, p_bool);

  // derived parameters from coarse file
  int gs = lastpt / save_pt;
  int num_steps = nsteps / save_step;
  double dr = (rmax - rmin) / ((double) lastpt);
  double dt = lam * dr;
  
  // bbhutil parameters & fields
  string fname0 = to_string(*resns[0]) + "-" + outfile;
  string fname1 = to_string(*resns[1]) + "-" + outfile;
  string fname2 = to_string(*resns[2]) + "-" + outfile;
  string phiname0 = pre1 + fname0 + ".sdf";
  string phiname1 = pre1 + fname1 + ".sdf";
  string phiname2 = pre1 + fname2 + ".sdf";
  string piname0 = pre2 + fname0 + ".sdf";
  string piname1 = pre2 + fname1 + ".sdf";
  string piname2 = pre2 + fname2 + ".sdf";
  bool wr2 = (pre2 == "0") ? false : true;
  int nwr = (wr2) ? 6 : 3; 
  char *name_arr[nwr];
  name_arr[0] = &phiname0[0], name_arr[1] = &phiname1[0],
    name_arr[2] = &phiname2[0];
  if (wr2) { name_arr[3] = &piname0[0], name_arr[4] = &piname1[0],
      name_arr[5] = &piname2[0]; }

  int npts0 = gs + 1;
  int npts1 = (same_grids) ? npts0 : 2*gs + 1;
  int npts2 = (same_grids) ? npts0 : 4*gs + 1;
  vector<double> phi_4h(npts0, 0.0), pi_4h((wr2) ? npts0 : 1, 0.0);
  vector<double> phi_2h(npts1, 0.0), pi_2h((wr2) ? npts1 : 1, 0.0);
  vector<double> phi_h(npts2, 0.0), pi_h((wr2) ? npts2 : 1, 0.0);
  double *field_arr[nwr];
  field_arr[0] = &phi_4h[0], field_arr[1] = &phi_2h[0],
    field_arr[2] = &phi_h[0];
  if (wr2) { field_arr[3] = &pi_4h[0], field_arr[4] = &pi_2h[0],
      field_arr[5] = &pi_h[0]; }
  vector<double> zeros(4, 0.0);
  vector<double> norms(4, 0.0);

  // output file
  ofstream ofs;
  ofs.open("conv-" + fname0 + ".csv", ofstream::out);
  ofs <<  "coarse,mid,fine,constant,points,times,same_times,same_grids\n"
      <<  phiname0+","+phiname1+","+phiname2+","+hold_const+"," << npts0
      <<","<< num_steps <<","<< boolalpha << same_times <<","<< same_grids
      << "\n\ndspn,dspn_bound,zero_pi,sommerfeld,tol,maxit\n" << dspn <<","
      << dspn_bound <<","<< zero_pi <<","<< sommerfeld <<","<< tol <<","
      << maxit << "\n\nr2m,rmin,rmax,ic_Dsq,ic_r0,ic_Amp\n" << r2m <<","
      << rmin <<","<< rmax <<","<< ic_Dsq <<","<< ic_r0 <<","<< ic_Amp
      << "\n\ncoarse grid:\nlastpt,save_pt,nsteps,save_step,\n" << lastpt
      <<","<< save_pt <<","<< nsteps <<","<< save_step << "\n\nlam,dr,dt,"
      << "tmax\n" << lam <<","<< dr <<","<< dt <<","<< dt*nsteps
      << "\n\ntime,Q" << pre1 << "(t),Q" << pre2 << "(t),sum" << endl;

  // iterate through time steps
  int t1 = ((same_times) ? 1 : 2);
  int t2 = ((same_times) ? 1 : 4);
  int r1 = ((same_grids) ? 1 : 2);
  int r2 = ((same_grids) ? 1 : 4);

  gft_set_multi();
  for (int t = 0; t < num_steps; ++t) {
    // compute ||u_4h(x,t) - u_2h(x,t)||/||u_2h(x,t) - u_h(x,t)||
    int times[nwr];
    times[0] = t+1, times[1] = t1*t+1, times[2] = t2*t+1;
    if (wr2) { times[3] = t+1, times[4] = t1*t+1, times[5] = t2*t+1; }
    read_step(name_arr, times, field_arr, nwr);
    norms = zeros;
    // iterate through grid points
    for (int j = 0; j < npts0; ++j) {
      norms[0] += abs(phi_4h[j] - phi_2h[r1*j]);
      norms[1] += abs(phi_2h[r1*j] - phi_h[r2*j]);
      if (wr2) {
	norms[2] += abs(pi_4h[j] - pi_2h[r1*j]);
	norms[3] += abs(pi_2h[r1*j] - pi_h[r2*j]);
      }
    }
    // write time, q_phi, q_pi, q_phi+pi
    ofs << t*save_step*dt <<","<< norms[0]/norms[1];
    if (wr2) { ofs <<","<< norms[2]/norms[3] <<","
		   << (norms[0] + norms[2])/(norms[1] + norms[3]); }
    ofs << endl;
  }
  gft_close_all();
  
  cout << "\nconv-" + fname0 + ".csv" << "  written with:" << endl;
  cout << "grid points used = " << npts0 << "  " << ((same_grids) ?
   "(same grids)" : "(dif grids)") << "\ntime steps used = " << num_steps
       << "  " << ((same_times) ? "(same times)" : "(dif times)") << endl;
  ofs.close();

  return 0;
}
