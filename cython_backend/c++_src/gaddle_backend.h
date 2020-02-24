#ifndef GADDLE_BACKEND_H
#define GADDLE_BACKEND_H

#include "gaddle_backend_global.h"
#include "common_libraries.h"

typedef std::vector<std::vector<double>> std_matrix;
std_matrix sqeuclidean_distance(std_matrix mat1, std_matrix mat2);

arma::mat convert_to_arma(std_matrix matrix);
std_matrix convert_to_std(arma::mat matrix);

std_matrix minimize_molecules(std_matrix mol1_pos, std_matrix mol2_pos,
                              std::vector<double> mol2_com, double sigma_scale,
                              long Nsteps, restraint_list restrains,
                              bonds bonds_info, bool same_com, double anchura,
                              std::vector<int> sim_types);

#endif // GADDLE_BACKEND_H
