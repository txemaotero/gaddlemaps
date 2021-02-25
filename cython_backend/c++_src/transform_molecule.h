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
#ifndef TRANSFORM_MOLECULE_H
#define TRANSFORM_MOLECULE_H
#include "common_libraries.h"

void move_mol_atom(arma::mat &atom_pos, bonds &bonds_info,
                   double sigma_scale=0.5);
arma::rowvec find_atom_random_displ(const arma::mat &atom_pos,
                                    bonds &bonds_info,
                                    int atom_index, double sigma_scale=0.5);
#endif // TRANSFORM_MOLECULE_H
