#ifndef ALGORITHM_UTILS_H
#define ALGORITHM_UTILS_H
#include "common_libraries.h"
#include <chrono>

class Chi2Calculator
{
public:
    Chi2Calculator(const arma::mat &molecule1, const arma::mat &molecule2,
                   const restraint_list &restrictions);
    double operator()(const arma::mat &new_mol2);
private:
    enum meth_to_call_enum {no_restrains, with_restrains, only_restrains};

    double chi2_molecules_restrains_contrib(const arma::mat &new_mol2);
    double chi2_molecules_only_restrains(const arma::mat &new_mol2);
    double chi2_molecules_with_restrains(const arma::mat &new_mol2);
    double chi2_molecules(const arma::mat &new_mol2);
    arma::mat mol1_positions;
    meth_to_call_enum meth_to_call;
    restraint_list restrictions;
    arma::uvec restriction2;
    int len_mol2;
    arma::mat mol1_not_restrictions;
    arma::mat mol1_restrictions;
};

bool accept_metropolis(double E0, double E1, double A=0.01);

arma::mat calc_sqeuclidean_distance(const arma::mat &matrix1,
                                    const arma::mat &matrix2);

class Progressbar
{
public:
    Progressbar(long total_steps, double start_chi);
    void operator()(long current_step, double current_chi);

private:
    long total_steps;
    long percentaje;
    long steps_done;
    double chi;
    void print();
    void end_print();
    std::chrono::steady_clock::time_point start;
};


class Progressbar_double
{
public:
    Progressbar_double(double total_steps);
    void operator()(double step_weight = 1);

private:
    double total_steps;
    long percentaje;
    double steps_done;
    void print();
    void end_print();
    std::chrono::steady_clock::time_point start;
};

#endif // ALGORITHM_UTILS_H
