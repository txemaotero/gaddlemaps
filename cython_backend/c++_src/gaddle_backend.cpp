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
#include "gaddle_backend.h"
#include "algorithm_utils.h"
#include "transform_molecule.h"

arma::mat convert_to_arma(std_matrix matrix){
    unsigned int long long nrows = matrix.size();
    unsigned int long long ncolumns = matrix[0].size();
    arma::mat matrixout(nrows, ncolumns);

    for (int i=0; i < nrows; ++i){
        for (int j=0; j < ncolumns; ++j){
            matrixout(i,j) = matrix[i][j];
        }
    }
    return matrixout;
}

std_matrix convert_to_std(arma::mat matrix){
    unsigned int long long nrows = matrix.n_rows;
    unsigned int long long ncolumns = matrix.n_cols;
    auto matrixout = std_matrix(nrows, std::vector<double>(ncolumns, 0.));
    for (int i=0; i < nrows; ++i){
        for (int j=0; j < ncolumns; ++j){
            matrixout[i][j] = matrix(i,j);
        }
    }
    return matrixout;
}


std_matrix sqeuclidean_distance(std_matrix mat1, std_matrix mat2){
    auto matrix1 = convert_to_arma(mat1);
    auto matrix2 = convert_to_arma(mat2);

    auto matrixout = calc_sqeuclidean_distance(matrix1, matrix2);
    return convert_to_std(matrixout);
}

std_matrix minimize_molecules(std_matrix mol1_pos, std_matrix mol2_pos,
                              std::vector<double> mol2_com, double sigma_scale,
                              long Nsteps, restraint_list restrains,
                              bonds bonds_info, bool same_com, double anchura,
                              std::vector<int> sim_types){

    arma::mat mol1 = convert_to_arma(mol1_pos);
    arma::mat mol2 = convert_to_arma(mol2_pos);

    arma::rowvec mol2com = arma::conv_to<arma::rowvec>::from(mol2_com);
    arma::mat out = minimize_molecules_routine(mol1, mol2, mol2com, sigma_scale, Nsteps,
                                      restrains, bonds_info, same_com, anchura,
                                      sim_types);
    mol2_pos = convert_to_std(out);
    return mol2_pos;
}
