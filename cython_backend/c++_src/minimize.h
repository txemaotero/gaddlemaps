#ifndef MINIMIZE_H
#define MINIMIZE_H
#include "common_libraries.h"

arma::mat minimize_molecules_routine(arma::mat &mol1_pos, arma::mat &mol2_pos,
                                arma::rowvec mol2_com, double sigma_scale,
                                long Nsteps, restraint_list &restrains,
                                bonds &bonds_info, bool same_com, double anchura,
                                const std::vector<int> sim_types);
#endif // MINIMIZE_H
