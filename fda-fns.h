#include <vector> // for everything
#include <cmath> // for ICs

using namespace std;

// x^2 and 1/x^2
inline double sq(double x) { return (x * x); }
inline double p4(double x) { return sq(sq(x)); }
inline double sqin(double x) { return (1.0 / sq(x)); }
inline double p4in(double x) { return sqin(sq(x)); }

// kreiss-oliger dissipation (p.23 choptuik notes)
inline double dissipate(double eps, const vector<double>& u, int ind)
{ return -eps * 0.0625 * ( u[ind-2] - 4*u[ind-1] + 6*u[ind]
			   - 4*u[ind+1] + u[ind+2] ); }

inline double symdiss1(double eps, const vector<double>& u)
{ return -eps * 0.0625 * ( u[3] - 4*u[2] + 7*u[1] - 4*u[0] ); }

inline double antidiss1(double eps, const vector<double>& u)
{ return -eps * 0.0625 * ( u[3] - 4*u[2] + 5*u[1] - 4*u[0] ); }

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
