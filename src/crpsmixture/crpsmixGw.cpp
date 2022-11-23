#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <cmath> 
#include <math.h>
#include <vector>

namespace py = pybind11;

double normpdf(double x) {
    return (1 / std::sqrt(2*M_PI)) * exp(-0.5 * x * x);
}

double auxcrps(double um, double h) {
    if (h == 0){
        return std::abs(um);
    } else {
        return 2 * h * normpdf(um / h) + um * std::erf(um / (h * std::sqrt(2))); 
    }
}

double crpsmixGw(py::array_t<double> m,py::array_t<double> w ,double y, double h) {
    py::buffer_info bufm = m.request();
    double *ptrm = static_cast<double *>(bufm.ptr);
    /*std::vector<double> array_m(m.size());
    std::memcpy(array_m.data(), m.data(), m.size()*sizeof(double));*/
    py::buffer_info bufw = w.request();
    double *ptrw = static_cast<double *>(bufw.ptr);

    
    py::ssize_t N = bufm.size;
   
    double crps1 = 0;
    double crps2 = 0;
    double W = 0;

    for (int i = 0; i < N; i++){
        W += ptrw[i];
	crps1 += ptrw[i] * auxcrps(y-ptrm[i], h);
        double crps3 = 0.5 * ptrw[i] * auxcrps(0.0, std::sqrt(2)*h);
	for (int j = 0; j < i; j++){
	    crps3 += ptrw[j] * auxcrps(ptrm[i] - ptrm[j], std::sqrt(2)*h);
	}
	crps2 += ptrw[i] * crps3;
    }
    return (crps1 - crps2 / W) / W;
}


PYBIND11_MODULE(_crpsmixture, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("crpsmixGw", &crpsmixGw, "Computes CRPS for mixture normal distribution for a one-dimensional observation y");
}
/*
<%
setup_pybind11(cfg)
%>
*/

