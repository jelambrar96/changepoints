#include "cost.h"



// ----------------------------------------------------------------------

MeanCost::MeanCost() {
    _n = 0;
}

MeanCost::~MeanCost() {
    if (_x_cusum != nullptr) {
        delete [] _x_cusum;
    }
    if (_x2_cusum != nullptr) {
        delete [] _x2_cusum;
    }
}

void MeanCost::setData(double *x, int n) {
    if (_n != n) {
        if (_x_cusum != nullptr) {
            delete [] _x_cusum;
        }
        if (_x2_cusum != nullptr) {
            delete [] _x2_cusum;
        }        
    }
    _n = n;
    if (_x_cusum == nullptr) {
        _x_cusum = new double[_n];
    }
    if (_x2_cusum == nullptr) {
        _x2_cusum = new double[_n];
    }

    double sum = 0, x2_sum = 0;
    for (int i = 0; i < _n; ++i) {
        sum += x[i];
        _x_cusum[i] = sum;
        x2_sum += (x[i] * x[i]);
        _x2_cusum[i] = x2_sum;
    }
}

double MeanCost::cost(int a, int b) {
    int m = b - a;
    double d_cusum = (_x_cusum[b] - _x_cusum[a]);
    double d2_cusum = (_x2_cusum[b] - _x2_cusum[a]);
    return d_cusum / m - d2_cusum / (m * m);
}

// ----------------------------------------------------------------------

LinearCost::LinearCost() {
    _n = 0;
}

LinearCost::~LinearCost() {
    if (_x_cusum != nullptr) {
        delete [] _x_cusum;
    }
    if (_x2_cusum != nullptr) {
        delete [] _x2_cusum;
    }
    if (_y_cusum != nullptr) {
        delete [] _y_cusum;
    }
    if (_y2_cusum != nullptr) {
        delete [] _y2_cusum;
    }
}

void LinearCost::setData(double *x, int n) {
    if (_n != n) {
        if (_x_cusum != nullptr) {
            delete [] _x_cusum;
        }
        if (_x2_cusum != nullptr) {
            delete [] _x2_cusum;
        }        
        if (_y_cusum != nullptr) {
            delete [] _y_cusum;
        }
        if (_y2_cusum != nullptr) {
            delete [] _y2_cusum;
        }        
    }
    _n = n;
    if (_x_cusum == nullptr) {
        _x_cusum = new double[_n];
    }
    if (_x2_cusum == nullptr) {
        _x2_cusum = new double[_n];
    }
    if (_y_cusum == nullptr) {
        _y_cusum = new double[_n];
    }
    if (_y2_cusum == nullptr) {
        _y2_cusum = new double[_n];
    }

    double x_sum = 0, x2_sum = 0;
    double y_sum = 0, y2_sum = 0;
    double xy_sum = 0;
    for (int i = 0; i < _n; ++i) {
        x_sum += (i+1);
        _x_cusum[i] = x_sum;
        x2_sum += ((i+1)*(i+1));
        _x2_cusum[i] = x2_sum;
        y_sum += (x[i]);
        _y_cusum[i] = y_sum;
        y2_sum += (x[i]*x[i]);
        _y2_cusum[i] = y2_sum;
        xy_sum = 0;
        _yx_cusum[i] = xy_sum;
    }
}

double LinearCost::cost(int a, int b) {
    int m = b - a;

    double y2cusum = _y2_cusum[b] - _y2_cusum[a];
    double ycusum = _y_cusum[b] - _y_cusum[a];

    double x2cusum = _x2_cusum[b] - _x2_cusum[a];
    double xcusum = _x_cusum[b] - _x_cusum[a];

    double xycusum = _yx_cusum[b] - _yx_cusum[a];

    double Syy = y2cusum - ycusum * ycusum / m;
    double Sxx = x2cusum - xcusum * xcusum / m;
    double Sxy = xycusum - xcusum * ycusum / m;

    return (Sxx == 0 ? std::numeric_limits<double>::max() : Syy - Sxy * Sxy / Sxx);
}

// ----------------------------------------------------------------------

StdCost::StdCost() {
    _n = 0;
}

StdCost::~StdCost() {
    if (_x_cusum != nullptr) {
        delete [] _x_cusum;
    }
    if (_x2_cusum != nullptr) {
        delete [] _x2_cusum;
    }
}

void StdCost::setData(double *x, int n) {
    if (_n != n) {
        if (_x_cusum != nullptr) {
            delete [] _x_cusum;
        }
        if (_x2_cusum != nullptr) {
            delete [] _x2_cusum;
        }        
    }
    _n = n;
    if (_x_cusum == nullptr) {
        _x_cusum = new double[_n];
    }
    if (_x2_cusum == nullptr) {
        _x2_cusum = new double[_n];
    }

    double sum = 0, x2_sum = 0;
    for (int i = 0; i < _n; ++i) {
        sum += x[i];
        _x_cusum[i] = sum;
        x2_sum += (x[i] * x[i]);
        _x2_cusum[i] = x2_sum;
    }
}

double StdCost::cost(int a, int b) {
    int m = b - a;
    double d_cusum = (_x_cusum[b] - _x_cusum[a]);
    double d2_cusum = (_x2_cusum[b] - _x2_cusum[a]);
    double ss = d2_cusum / m - d_cusum * d_cusum / (m * m);
    double c = (ss > 0) ? ss : std::numeric_limits<double>::epsilon();
    return m * log(c);
}

// ----------------------------------------------------------------------

RmsCost::RmsCost() {
    _n = 0;
}

RmsCost::~RmsCost() {
    if (_x2_cusum != nullptr) {
        delete [] _x2_cusum;
    }
}

void RmsCost::setData(double *x, int n) {
    if (_n != n) {
        if (_x2_cusum != nullptr) {
            delete [] _x2_cusum;
        }        
    }
    _n = n;
    if (_x2_cusum == nullptr) {
        _x2_cusum = new double[_n];
    }

    double x2_sum = 0;
    for (int i = 0; i < _n; ++i) {
        x2_sum += (x[i] * x[i]);
        _x2_cusum[i] = x2_sum;
    }
}

double RmsCost::cost(int a, int b) {
    int m = b - a;
    double d2_cusum = (_x2_cusum[b] - _x2_cusum[a]);
    double c = (d2_cusum > 0) ? d2_cusum : std::numeric_limits<double>::epsilon();
    return m * (log(c) - log(m));
}

// ----------------------------------------------------------------------


