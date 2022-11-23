#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath> 
#include <math.h>
#include <vector>


#include "thirdparty/boost/boost/math/quadrature/gauss_kronrod.hpp"
#include "thirdparty/boost/boost/math/distributions/students_t.hpp"


#include <iostream>
using boost::math::quadrature::gauss_kronrod;
using boost::math::students_t;

namespace py = pybind11;

double hypergeometric( double a, double b, double c, double x )
{
   const double TOLERANCE = 1.0e-10;
   double term = a * b * x / c;
   double value = 1.0 + term;
   int n = 1;

   while ( abs( term ) > TOLERANCE )
   {
      a++, b++, c++, n++;
      term *= a * b * x / c / n;
      value += term;
   }

   return value;
}



double tcdf0(double x, double df, double m, double h) {
    double c = 3 / 2;
    double a = 0.5;
    double b = 0.5*(df + 1);
    double z = -1*(x*x) / df;
    double HG = hypergeometric(a, b, c, z);
    double T = HG / (std::sqrt(df * M_PI) * std::tgamma(0.5*df));
    return (0.5 + x * std::tgamma(b)* T);
}

double tcdf(double x, double df, double m, double h){
    students_t dist(df);
    double z = (x - m) / h;
    return cdf(dist, z);
}


double cpp_int_lims(py::array_t<double>& Y, py::array_t<double>& M, py::array_t<double>& W, py::array_t<int>& I, double h, double df, double low, double up) {

    py::buffer_info bufY = Y.request();
    double* ptrY = (double*)bufY.ptr;
    py::ssize_t n = bufY.shape[0];

    py::buffer_info bufM = M.request();
    double* ptrM = (double*)bufM.ptr;
    
    py::buffer_info bufW = W.request();
    double* ptrW = (double*)bufW.ptr;

    py::buffer_info bufI = I.request();
    int* ptrI = (int*)bufI.ptr;
 
    double out = 0;
        
    double inf = std::numeric_limits<double>::infinity();
    
    for (int i = 0; i < n; i++){
        double yval = ptrY[i];
        int s0 = ptrI[i];
        int s1 = ptrI[i+1];
        int end = s1 - s0;
        std::vector<double> w0(end);
        std::vector<double> m0(end);
        for(int k = 0; k < end; k++){
            w0[k] = ptrW[k+s0];
            m0[k] = ptrM[k+s0];
        }        
   
        students_t dist(df);
       
        auto f1 = [&dist, &h, &m0, &w0, &end](double t) 
            { 
		double sum = 0;
 		for(int j = 0; j < end; j++){
		    double z = (t-m0[j])/h;
		    sum = sum + w0[j]*cdf(dist, z);
		}
                return sum*sum;		
            };
        auto f2 = [&dist, &h, &m0, &w0, &end](double t) 
            { 
	        double sum = 0;
 		for(int j = 0; j < end; j++){
		    double z = (t-m0[j])/h;
		    sum = sum + w0[j]*cdf(dist, z);
		}

		return (1-sum)*(1-sum); 
	    };
        double error;
        double I1 = gauss_kronrod<double, 15>::integrate(f1, low, yval, 0, 0, &error);
        double I2 = gauss_kronrod<double, 15>::integrate(f2, yval, up, 0, 0, &error);

        out = out + I1 + I2; 
    
    }
    return out / n;
}


PYBIND11_MODULE(_crpsmixture, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("cpp_int_lims", &cpp_int_lims, "Computes CRPS for mixture normal distribution for a one-dimensional observation y");
}
/*
<%
setup_pybind11(cfg)
%>
*/

