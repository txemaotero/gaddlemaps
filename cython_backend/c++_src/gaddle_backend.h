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
