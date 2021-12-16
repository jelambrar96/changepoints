#include "cost.h"



// ----------------------------------------------------------------------

MeanCost::MeanCost() {
    _n = 0;
    _x_cusum = nullptr;
    _x2_cusum = nullptr;
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
    double d_cusum = (a == -1) ? _x_cusum[b] : (_x_cusum[b] - _x_cusum[a]);
    double d2_cusum = (a == -1) ? _x2_cusum[b] : (_x2_cusum[b] - _x2_cusum[a]);
    return d_cusum / m - d2_cusum / (m * m);
}

double MeanCost::getLmin() { return 1; }

double MeanCost::getFWD(double * d) {
    for (int i = 0; i < _n; ++i) {
        d[i] = _x2_cusum[i] - _x_cusum[i] * _x_cusum[i] / (i+1);
    }
}

double MeanCost::getREV(double * rev) {
    int n = _n;
    for (int i = n - 1; i >= 0; --i) {
        double cusum = _x_cusum[n-1] - _x_cusum[i];
        double cusum2 = _x2_cusum[n-1] - _x2_cusum[i];
        int temp_n = n - i;
        rev[i] = cusum2 - cusum * cusum / temp_n;
    }
}

// ----------------------------------------------------------------------

LinearCost::LinearCost() {
    _n = 0;
    _x_cusum = nullptr;
    _x2_cusum = nullptr;
    _y_cusum = nullptr;
    _y2_cusum = nullptr;
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

    double y2cusum = (a==-1) ? _y2_cusum[b] : _y2_cusum[b] - _y2_cusum[a];
    double ycusum = (a==-1) ? _y_cusum[b] : _y_cusum[b] - _y_cusum[a];

    double x2cusum = (a==-1) ? _x2_cusum[b] : _x2_cusum[b] - _x2_cusum[a];
    double xcusum =  (a==-1) ? _x_cusum[b] : _x_cusum[b] - _x_cusum[a];

    double xycusum = (a==-1) ? _yx_cusum[b] : _yx_cusum[b] - _yx_cusum[a];

    double Syy = y2cusum - ycusum * ycusum / m;
    double Sxx = x2cusum - xcusum * xcusum / m;
    double Sxy = xycusum - xcusum * ycusum / m;

    return (Sxx == 0 ? std::numeric_limits<double>::max() : Syy - Sxy * Sxy / Sxx);
}

double LinearCost::getLmin() { return 2; }

double LinearCost::getFWD(double * fwd) {
    int n = _n;
    for (int i = 0; i < n; ++i) {
        double ycusum = _y_cusum[i];
        double xcusum = _x_cusum[i];
        double x2cusum = _x2_cusum[i];
        double y2cusum = _y2_cusum[i];
        double xycusum = _yx_cusum[i];
        double Syy = y2cusum - ycusum * ycusum / (i + 1);
        double Sxx = x2cusum - xcusum * xcusum / (i + 1);
        double Sxy = xycusum - xcusum * ycusum / (i + 1);
        fwd[i] = (Sxx == 0 ? std::numeric_limits<double>::max() : Syy - Sxy * Sxy / Sxx);
    }    
}

double LinearCost::getREV(double * rev) { 
    int n = _n;
    for (int i = n - 1; i >= 0; --i) {
        int temp_n = n - i;
        double ycusum = _y_cusum[n-1] - _y_cusum[i];
        double xcusum = _x_cusum[n-1] - _x_cusum[i];
        double x2cusum = _x2_cusum[n-1] - _x2_cusum[i];
        double y2cusum = _y2_cusum[n-1] - _y2_cusum[i];
        double xycusum = _yx_cusum[n-1] - _yx_cusum[i];
        double Syy = y2cusum - ycusum * ycusum / (n - i);
        double Sxx = x2cusum - xcusum * xcusum / (n - i);
        double Sxy = xycusum - xcusum * ycusum / (n - i);
        rev[i] = (Sxx == 0 ? std::numeric_limits<double>::max() : Syy - Sxy * Sxy / Sxx);
    }
}

// ----------------------------------------------------------------------

StdCost::StdCost() {
    _n = 0;
    _x_cusum = nullptr;
    _x2_cusum = nullptr;
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
    double d_cusum = (a == -1) ? _x_cusum[b] : (_x_cusum[b] - _x_cusum[a]);
    double d2_cusum = (a == -1) ? _x2_cusum[b] : (_x2_cusum[b] - _x2_cusum[a]);
    double ss = d2_cusum / m - d_cusum * d_cusum / (m * m);
    double c = (ss > 0) ? ss : std::numeric_limits<double>::epsilon();
    return m * log(c);
}

double StdCost::getLmin() { return 2; }

double StdCost::getFWD(double * fwd) {
    int n = _n;
    for (int i = 0; i < n; ++i) {
        double cusum = _x_cusum[i];
        double cusum2 = _x2_cusum[i];
        int temp_n = i + 1;
        double ss = cusum2 / temp_n - cusum * cusum / (temp_n * temp_n);
        double c = (ss > 0) ? ss : std::numeric_limits<double>::epsilon(); 
        sumlog = log(c);
        fwd[i] = temp_n * sumlog;
    }
}

double StdCost::getREV(double * rev) {
    int n = _n;
    for (int i = n - 1; i >= 0; --i) {
        double cusum = _x_cusum[n-1] - _x_cusum[i];
        double cusum2 = _x2_cusum[n-1] - _x2_cusum[i];
        int temp_n = n - i;
        double ss = cusum2 / temp_n - cusum * cusum / (temp_n * temp_n);
        double c = (ss > 0) ? ss : std::numeric_limits<double>::epsilon(); 
        sumlog = log(c);
        rev[i] = temp_n * sumlog;
    }
}

// ----------------------------------------------------------------------

RmsCost::RmsCost() {
    _n = 0;
    _x_cusum = nullptr;
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
    double d2_cusum = (a==-1) ? _x2_cusum[b] : (_x2_cusum[b] - _x2_cusum[a]);
    double c = (d2_cusum > 0) ? d2_cusum : std::numeric_limits<double>::epsilon();
    return m * (log(c) - log(m));
}

double RmsCost::getLmin() { return 2; }

double RmsCost::getFWD(double * fwd) {
    int n = _n;
    for (int i = 0; i < n; ++i) {
        // cusum
        int temp_n = i + 1;
        double cusum2 = _x2_cusum[i];
        double c = (cusum2 > 0) ? cusum2 : std::numeric_limits<double>::epsilon();
        sumlog = (log(c) - log(temp_n));
        fwd[i] = temp_n * sumlog;
    }
}

double RmsCost::getREV(double * d) {
    int n = _n;
    for (int i = n - 1; i >= 0; --i) {
        int temp_n = n - i;
        double cusum2 = _x2_cusum[n-1] - _x2_cusum[i];
        double c = (cusum2 > 0) ? cusum2 : std::numeric_limits<double>::epsilon(); 
        sumlog = (log(c) - log(temp_n));
        rev[i] = temp_n * sumlog;
    }
}

// ----------------------------------------------------------------------


