#ifndef RESIDUALS_H
#define RESIDUALS_H

#include <cmath>
#include <limits>

// #include <iostream>


int getResidueIndex(double * fwd, double * rev, int n, int Lmin);
int argmin(double * d, int n);

void linearResidual(double *y, double * fwd, double * rev, int n);
void linearResidual(double *y, double * fwd, int n);

void meanResidual(double * x, double * fwd, double * rev, int n);
void meanResidual(double * x, double * fwd, int n);

void rmsResidual(double * x, double * fwd, double * rev, int n);
void rmsResidual(double * x, double * fwd, int n);

void stdResidual(double * x, double * fwd, double * rev, int n);
void stdResidual(double * x, double * fwd, int n);


#endif

