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
#ifndef MINIMIZE_H
#define MINIMIZE_H
#include "common_libraries.h"

arma::mat minimize_molecules_routine(arma::mat &mol1_pos, arma::mat &mol2_pos,
                                arma::rowvec mol2_com, double sigma_scale,
                                long Nsteps, restraint_list &restrains,
                                bonds &bonds_info, bool same_com, double anchura,
                                const std::vector<int> sim_types);
#endif // MINIMIZE_H
