#ifndef TRANSFORM_MOLECULE_H
#define TRANSFORM_MOLECULE_H
#include "common_libraries.h"

void move_mol_atom(arma::mat &atom_pos, bonds &bonds_info,
                   double sigma_scale=0.5);
arma::rowvec find_atom_random_displ(const arma::mat &atom_pos,
                                    bonds &bonds_info,
                                    int atom_index, double sigma_scale=0.5);
#endif // TRANSFORM_MOLECULE_H
