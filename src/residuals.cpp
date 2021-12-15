#include "residuals.h"


int getResidueIndex(double * fwd, double * rev, int n, int Lmin) {

    int m = n - 2 * Lmin;
    double * x = new double[m];
    int LminN = n - Lmin;
    for (int i = 0; i < m; ++i) {
        x[i] = fwd[i + Lmin] + rev[i + Lmin];
    }
    int xargmin =  argmin(x, m);
    delete [] x;
    return xargmin + Lmin;
}

int argmin(double * d, int n) {
    double min = d[0];
    double argmin = 0;
    for (int i = 0; i < n; ++i) {
        if (d[i] < min) {
            min = d[i];
            argmin = i;
        }
    }
    return argmin;
}


void linearResidual(double *y, double * fwd, double * rev, int n) {

    double xcusum = 0;
    double ycusum = 0;
    double x2cusum = 0;
    double y2cusum = 0;
    double xycusum = 0;

    for (int i = 0; i < n; ++i) {
        
        ycusum += y[i];
        xcusum += (i + 1);

        x2cusum += ((i + 1) * (i + 1));
        y2cusum += (y[i] * y[i]);
        xycusum += ((i + 1) * y[i]);

        double Syy = y2cusum - ycusum * ycusum / (i + 1);
        double Sxx = x2cusum - xcusum * xcusum / (i + 1);
        double Sxy = xycusum - xcusum * ycusum / (i + 1);

        fwd[i] = (Sxx == 0 ? std::numeric_limits<double>::max() : Syy - Sxy * Sxy / Sxx);
    }

    xcusum = 0;
    ycusum = 0;
    x2cusum = 0;
    y2cusum = 0;
    xycusum = 0;

    for (int i = n - 1; i >= 0; --i) {

        int temp_n = n - i;

        ycusum += y[i];
        xcusum += (i + 1);

        x2cusum += ((i + 1) * (i + 1));
        y2cusum += (y[i] * y[i]);
        xycusum += ((i + 1) * y[i]);

        double Syy = y2cusum - ycusum * ycusum / (n - i);
        double Sxx = x2cusum - xcusum * xcusum / (n - i);
        double Sxy = xycusum - xcusum * ycusum / (n - i);

        rev[i] = (Sxx == 0 ? std::numeric_limits<double>::max() : Syy - Sxy * Sxy / Sxx);

    }

}



void linearResidual(double *y, double * fwd, int n) {

    double xcusum = 0;
    double ycusum = 0;
    double x2cusum = 0;
    double y2cusum = 0;
    double xycusum = 0;

    for (int i = 0; i < n; ++i) {
        
        ycusum += y[i];
        xcusum += (i + 1);

        x2cusum += ((i + 1) * (i + 1));
        y2cusum += (y[i] * y[i]);
        xycusum += ((i + 1) * y[i]);

        double Syy = y2cusum - ycusum * ycusum / (i + 1);
        double Sxx = x2cusum - xcusum * xcusum / (i + 1);
        double Sxy = xycusum - xcusum * ycusum / (i + 1);

        fwd[i] = Syy;
        fwd[i] += (Sxx == 0 ? std::numeric_limits<double>::max() : Sxy * Sxy / Sxx);
    }

}


void meanResidual(double * x, double * fwd, double * rev, int n) {
    
    // auxiliar
    double cusum;
    double cusum2;

    // compute cusum
    cusum = 0;
    cusum2 = 0;

    for (int i = 0; i < n; ++i) {
        // cusum
        cusum += x[i];
        cusum2 += x[i] * x[i];
        fwd[i] = cusum2 - cusum * cusum / (i+1);
    }

    cusum = 0;
    cusum2 = 0;

    for (int i = n - 1; i >= 0; --i) {
        cusum += x[i];
        cusum2 += x[i] * x[i];
        int temp_n = n - i;
        rev[i] = cusum2 - cusum * cusum / temp_n;
    }

}

void meanResidual(double * x, double * fwd, int n) {
    
    // auxiliar
    double cusum;
    double cusum2;

    // compute cusum
    cusum = 0;
    cusum2 = 0;

    for (int i = 0; i < n; ++i) {
        // cusum
        cusum += x[i];
        cusum2 += x[i] * x[i];
        fwd[i] = cusum2 - cusum * cusum / (i+1);
    }

}


void rmsResidual(double * x, double * fwd, double * rev, int n) {

    // compute cusum
    double cusum2 = 0;
    double sumlog;
    // double logn = log(n);

    cusum2 = 0;
    sumlog = 0;
    for (int i = 0; i < n; ++i) {
        // cusum
        cusum2 += (x[i] * x[i]);
        double c = (cusum2 > 0) ? cusum2 : std::numeric_limits<double>::epsilon();
        sumlog = (log(c) - log(i + 1));
        fwd[i] = (i + 1) * sumlog;
    }

    cusum2 = 0;
    sumlog = 0;
    for (int i = n - 1; i >= 0; --i) {
        int temp_n = n - i;
        cusum2 += (x[i] * x[i]);
        double c = (cusum2 > 0) ? cusum2 : std::numeric_limits<double>::epsilon(); 
        sumlog = (log(c) - log(temp_n));
        rev[i] = temp_n * sumlog;
    }

}

void rmsResidual(double * x, double * fwd, int n) {

    // compute cusum
    double cusum2;
    double sumlog;
    // double logn = log(n);

    cusum2 = 0;
    sumlog = 0;
    for (int i = 0; i < n; ++i) {
        // cusum
        cusum2 += (x[i] * x[i]);
        double c = (cusum2 > 0) ? cusum2 : std::numeric_limits<double>::epsilon();
        sumlog = (log(c) - log(i + 1));
        fwd[i] = (i + 1) * sumlog;
    }

}

void stdResidual(double * x, double * fwd, double * rev, int n) {

    double cusum2; //  = 0;
    double cusum; //  = 0;
    double sumlog;

    cusum2 = 0;
    cusum = 0;
    sumlog = 0;
    for (int i = 0; i < n; ++i) {
        cusum2 += (x[i] * x[i]);
        cusum += x[i];
        double ss = cusum2 / (i + 1) - cusum * cusum / ((i + 1) * (i + 1));
        double c = (ss > 0) ? ss : std::numeric_limits<double>::epsilon(); 
        sumlog = log(c);
        fwd[i] = (i + 1) * sumlog;
    }

    cusum2 = 0;
    cusum = 0;
    sumlog = 0;
    for (int i = (int)n - 1; i >= 0; --i) {
        int temp_n = n - i;
        cusum2 += (x[i] * x[i]);
        cusum += x[i];
        double ss = cusum2 / temp_n - cusum * cusum / (temp_n * temp_n);
        double c = (ss > 0) ? ss : std::numeric_limits<double>::epsilon(); 
        sumlog = log(c);
        rev[i] = (temp_n) * sumlog;
    }

}


void stdResidual(double * x, double * fwd, int n) {

    double cusum2; //  = 0;
    double cusum; //  = 0;
    double sumlog;

    cusum2 = 0;
    cusum = 0;
    sumlog = 0;
    for (int i = 0; i < n; ++i) {
        cusum2 += (x[i] * x[i]);
        cusum += x[i];
        double ss = cusum2 / (i + 1) - cusum * cusum / ((i + 1) * (i + 1));
        double c = (ss > 0) ? ss : std::numeric_limits<double>::epsilon(); 
        sumlog = log(c);
        fwd[i] = (i + 1) * sumlog;
    }

}