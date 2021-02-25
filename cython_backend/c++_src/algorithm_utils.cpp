//    Gaddlemaps python module.
//    Copyright (C) 2019-2021 José Manuel Otero Mato, Hadrián Montes Campos, Luis Miguel Varela Cabo
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with/ this program.  If not, see <https://www.gnu.org/licenses/>.
#include "algorithm_utils.h"

auto double_dist = std::uniform_real_distribution<double>(0,1.);
auto random_engine1 = std::default_random_engine();


double chi2_molecules(const arma::mat &molecule1, const arma::mat &molecule2,
                      const restraint_list &restraints){
    // The original code checks that the restraints list is not empty
    // and that molecule2 is larger than molecule 1. In this implementation
    // We asume that the parameters were introduced in the correct order.

    auto dist_matrix = calc_sqeuclidean_distance(molecule1, molecule2);
    arma::colvec d_min_row = arma::min(dist_matrix, 1);

    // Update the values of the minimum distances acording to the constraints

    for (int restriction = 0; restriction < restraints.size(); ++restriction){
        auto rest = restraints[restriction];
        d_min_row(rest[0]) = dist_matrix(rest[0], rest[1]);
    }
    double chi2 = arma::accu(d_min_row);
    arma::ucolvec NCG_far_vec = arma::unique(arma::index_min(dist_matrix, 1));
    int NCG_far = int(molecule2.n_rows - NCG_far_vec.n_elem);

    // Multiplicative factor to enforce restraints
    double factor = 1;
    for (int i=0; i <= NCG_far; ++i){
        factor *= 1.1;
    }
    return chi2*factor;
}

Chi2Calculator::Chi2Calculator(const arma::mat &molecule1, const arma::mat &molecule2,
                               const restraint_list &restrictions)
{
    mol1_positions = molecule1;
    if (restrictions.empty())
    {
        std::cout << "empty\n";
        meth_to_call = meth_to_call_enum::no_restrains;
    }
    else{
        meth_to_call = meth_to_call_enum::with_restrains;

        this->restrictions = restrictions;
        auto restriction1 = std::vector<arma::u64>();
        auto restriction2_aux = std::vector<arma::u64>();
        for (auto &rest:restrictions)
        {
            restriction1.push_back(rest[0]);
            restriction2_aux.push_back(rest[1]);
        }

        len_mol2 = molecule2.n_rows;

        auto not_restriction1 = std::vector<arma::u64>();
        for (int i=0; i < molecule1.n_rows; ++i)
        {
            auto index = std::find(restriction1.begin(), restriction1.end(), i);
            if (index == restriction1.end())
            {
                not_restriction1.push_back(i);
            }
        }
        mol1_restrictions = molecule1.rows(arma::uvec(restriction1));
        mol1_not_restrictions = molecule1.rows(arma::uvec(not_restriction1));
        if (not_restriction1.empty())
        {
            meth_to_call = meth_to_call_enum::only_restrains;
            // std::set<int> restriction2_set(restriction2_aux.begin(), restriction2_aux.end());
            // n_cg_far_fact = std::pow(1.1, len_mol2 - restriction2_set.size());
        }
        restriction2 = arma::uvec(restriction2_aux);
    }
}

double Chi2Calculator::chi2_molecules_restrains_contrib(const arma::mat &new_mol2)
{
    arma::mat mol2_restrains = new_mol2.rows(restriction2);
    return arma::accu(arma::square(mol2_restrains - mol1_restrictions));
}

double Chi2Calculator::chi2_molecules_only_restrains(const arma::mat &new_mol2)
{
    return chi2_molecules_restrains_contrib(new_mol2);
}

double Chi2Calculator::chi2_molecules_with_restrains(const arma::mat &new_mol2)
{
    auto chi2 = chi2_molecules_restrains_contrib(new_mol2);

    auto dist_matrix = calc_sqeuclidean_distance(mol1_not_restrictions, new_mol2);

    arma::colvec d_min_row = arma::min(dist_matrix, 1);

    // Update the values of the minimum distances acording to the constraints

    chi2 += arma::accu(d_min_row);
    auto min_index = arma::index_min(dist_matrix, 1);
    arma::ucolvec NCG_far_vec = arma::unique(arma::join_cols(min_index, restriction2));

    int NCG_far = int(len_mol2 - NCG_far_vec.n_elem);

    // Multiplicative factor to enforce restraints
    double factor = 1;
    for (int i=0; i <= NCG_far; ++i){
        factor *= 1.1;
    }
    return chi2*factor;
}

double Chi2Calculator::chi2_molecules(const arma::mat &new_mol2)
{
    auto dist_matrix = calc_sqeuclidean_distance(mol1_positions, new_mol2);

    arma::colvec d_min_row = arma::min(dist_matrix, 1);

    // Update the values of the minimum distances acording to the constraints

    auto chi2 = arma::accu(d_min_row);
    arma::ucolvec NCG_far_vec = arma::unique(arma::index_min(dist_matrix, 1));

    int NCG_far = int(len_mol2 - NCG_far_vec.n_elem);

    // Multiplicative factor to enforce restraints
    double factor = 1;
    for (int i=0; i <= NCG_far; ++i){
        factor *= 1.1;
    }
    return chi2*factor;

}

double Chi2Calculator::operator()(const arma::mat &new_mol2)
{
    switch (meth_to_call)
    {
        case meth_to_call_enum::no_restrains:
            return chi2_molecules(new_mol2);

        case meth_to_call_enum::with_restrains:
            return chi2_molecules_with_restrains(new_mol2);
    
        case meth_to_call_enum::only_restrains:
            return chi2_molecules_only_restrains(new_mol2);
    }
}

bool accept_metropolis(double E0, double E1, double A)
{
    auto R = E0/E1;
    if (R >= 1){
        return true;
    }
    else{
        bool cond = (double_dist(random_engine1) <= A*R);
        return cond;
    }
}

arma::mat calc_sqeuclidean_distance(const arma::mat &matrix1,
                                    const arma::mat &matrix2){
    int N = int(matrix1.n_rows);
    int K = int(matrix2.n_rows);
    arma::mat XX(N,1);
    arma::mat YY(1,K);
    arma::mat XY (N,K);
    arma::mat D (N,K);

    XX = arma::sum((matrix1 % matrix1), 1);
    YY = arma::sum((matrix2 % matrix2), 1).t();
    XY = (2*matrix1)*matrix2.t();

    arma::mat ONES1 = arma::ones<arma::mat>(1, K);
    arma::mat ONES2 = arma::ones<arma::mat>(N, 1);

    D = XX*ONES1 + ONES2*YY - XY;

    return D;
}

Progressbar::Progressbar(long total_steps, double start_chi)
{
    this->total_steps = total_steps;
    percentaje = -1;
    chi = start_chi;
    steps_done = 0;
    start = std::chrono::steady_clock::now();
}

void Progressbar::operator ()(long current_step, double current_chi){
    steps_done = current_step;
    if (current_chi != chi)
    {
        chi = current_chi;
        print();
    }
    if (current_step == 0){
        start = std::chrono::steady_clock::now();
    }
    if ((100 * steps_done / total_steps) != percentaje){
        percentaje = 100 * steps_done / total_steps;
        print();
    }
    if (steps_done == total_steps){
        end_print();
    }
}

void Progressbar::print(){
    std::printf("%3i", int(percentaje));
    std::cout << "% eta:";
    std::chrono::steady_clock::time_point current =std::chrono::steady_clock::now();
    if (percentaje !=0){
        auto elapsed = std::chrono::duration_cast <std::chrono::microseconds>(current-start).count();
        elapsed = elapsed * (double(100./percentaje) - 1)/1000000;
        int seconds = elapsed % 60;
        int minutes = elapsed / 60 % 60;
        int hours = int(elapsed /3600);
        std::printf(" %5ih %2imin %2is", hours, minutes, seconds);
    }
    else{
        std::cout << "   ---h --min --s";
    }
    std::cout << "  step " << steps_done << "/" << total_steps;
    std::cout << " Chi2 = " << chi;
    std::cout << "\r";
    std::cout.flush();
}

void Progressbar::end_print(){
    std::printf("%3i", int(percentaje));
    std::cout << "% elapsed time:";
    std::chrono::steady_clock::time_point current = std::chrono::steady_clock::now();

    auto elapsed = std::chrono::duration_cast <std::chrono::microseconds>(current-start).count();
    elapsed = elapsed / 1000000;
    int seconds = elapsed % 60;
    int minutes = elapsed / 60 % 60;
    int hours = int(elapsed /3600);
    std::printf(" %5ih %2imin %2is", hours, minutes, seconds);
    std::cout << " Chi2 = " << chi;
    std::cout << std::endl;;
}

Progressbar_double::Progressbar_double(double total_steps)
{
    this->total_steps = total_steps;
    percentaje = -1;
    steps_done = 0;
    start = std::chrono::steady_clock::now();
}

void Progressbar_double::operator ()(double step_weight){
    steps_done += step_weight;

    if (long(100. * steps_done / total_steps) != percentaje){
        percentaje = 100. * steps_done / total_steps;
        print();
    }
    if (steps_done >= total_steps){
        end_print();
    }
}

void Progressbar_double::print(){
    std::printf("%3i", int(percentaje));
    std::cout << "% eta:";
    std::chrono::steady_clock::time_point current =std::chrono::steady_clock::now();
    if (percentaje !=0){
        auto elapsed = std::chrono::duration_cast <std::chrono::microseconds>(current-start).count();
        elapsed = elapsed * (double(100./percentaje) - 1)/1000000;
        int seconds = elapsed % 60;
        int minutes = elapsed / 60 % 60;
        int hours = int(elapsed /3600);
        std::printf(" %5ih %2imin %2is", hours, minutes, seconds);
    }
    else{
        std::cout << "   ---h --min --s";
    }
    std::cout << "\r";
    std::cout.flush();
}

void Progressbar_double::end_print(){
    std::printf("%3i", int(percentaje));
    std::cout << "% elapsed time:";
    std::chrono::steady_clock::time_point current = std::chrono::steady_clock::now();

    auto elapsed = std::chrono::duration_cast <std::chrono::microseconds>(current-start).count();
    elapsed = elapsed / 1000000;
    int seconds = elapsed % 60;
    int minutes = elapsed / 60 % 60;
    int hours = int(elapsed /3600);
    std::printf(" %5ih %2imin %2is", hours, minutes, seconds);

    std::cout << std::endl;;
}
