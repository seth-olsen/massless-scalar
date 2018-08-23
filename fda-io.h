#include <string> // for parameter input
#include <map> // for parameter input
#include <cstdlib> // for atoi() and atof()
#include "bbhutil.h" // for output to .sdf

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

// write fields using bbhutil
void write_step(double* fields[], int nfields, char* files[],
		double time, int* shape, int rank, double* coordinates) {
  for (int k = 0; k < nfields; ++k) {
    gft_out_bbox(files[k], time, shape, rank, coordinates, fields[k]);
  }
  return;
}



