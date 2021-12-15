#include "perl.h"


std::vector<int> perl(AbstractCost * cost; double pen = NULL) {

    int n = cost->getN();

    double penin = (pen == NULL) ? log(n) : pen;

    double * F = new double[n + 1];
    for (int i = 0; i < n +1; ++i) {
        F[i] = 0;
    }

    std::vector<double> R;
    R.push_back(0);

    int * candidates = new int[n + 1];
    for (int i = 0; i < n +1; ++i) {
        candidates[i] = 0;
    }

    F[0] -= penin;
    // int size_r = 1;

    for (int i = 2; i < n + 1; ++i) {

        int size_r = (int)R.size();
        double seg_costs = new double[size_r];

        for (int j = 0; j < size_r) {
            seg_costs[j] = cost->cost(R[j], i); 
        }

        double * Fcost = new double[size_r];
        for (int j = 0; j < size_r; ++j) {
            Fcost[j] = F[R[j]] + seg_costs[j];
        }

        int tau = argmin(Fcost, size_r);
        F[i] = Fcost[tau];

        candidates[i] = R[tau];

        R.resize(0); // free R
        for (int j=0; j < size_r; ++j) {
            if (F[j] < F[i]) {
                R.push_back(j);
            }
        }
        R.push_back(i - 1);

        delete [] Fcost;
        delete [] seg_costs;

    }

    std::vector<int> changepoints;
    int last = candidates[n];
    changepoints.push_back(last);
    while (last > 0) {
        last = candidates[last];
        changepoints.push_back(last);
    }

    delete [] F;
    delete [] candidates;

    std::sort(changepoints.begin(), changepoints.end());

    return changepoints;

}


