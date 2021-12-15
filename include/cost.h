#ifndef COST_H
#define COST_H


#include <cmath>
#include <limits>


class AbstractCost {
public:
    AbstractCost() {};
    virtual ~AbstractCost() {};
    virtual void setData(double *d, int n) = 0;
    virtual double cost(int a, int b) = 0;
    virtual int getN() { return _n; }
protected:
    int _n;
};


class MeanCost : public AbstractCost {
public:
    MeanCost();
    virtual ~MeanCost();
    virtual void setData(double *d, int n);
    virtual double cost(int a, int b);

private:
    double * _x_cusum;
    double * _x2_cusum;
};


class LinearCost : public AbstractCost {
public:
    LinearCost();
    virtual ~LinearCost();
    virtual void setData(double *d, int n);
    virtual double cost(int a, int b);

private:
    double * _x_cusum;
    double * _x2_cusum;
    double * _y_cusum;
    double * _y2_cusum;
    double * _yx_cusum;
};


class StdCost : public AbstractCost {
public:
    StdCost();
    virtual ~StdCost();
    virtual void setData(double *d, int n);
    virtual double cost(int a, int b);

private:
    double * _x_cusum;
    double * _x2_cusum;

};


class RmsCost : public AbstractCost {
public:
    RmsCost();
    virtual ~RmsCost();
    virtual void setData(double *d, int n);
    virtual double cost(int a, int b);

private:
    double * _x_cusum;
    double * _x2_cusum;
};

#endif
