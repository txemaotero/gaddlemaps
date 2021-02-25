#    Gaddlemaps python module.
#    Copyright (C) 2019-2021 José Manuel Otero Mato, Hadrián Montes Campos, Luis Miguel Varela Cabo
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.
def py_squeclidean_distance(matrix1, matrix2):
    return py_squeclidean_distance_back(matrix1, matrix2)

    
cdef std_matrix py_squeclidean_distance_back(std_matrix matrix1,
                                             std_matrix matrix2):
    return sqeuclidean_distance(matrix1, matrix2)

def py_minimize_molecules(mol1_pos, mol2_pos, mol2_com, sigma_scale, Nsteps,
                          restrains, bonds_info, anchura, sim_types):
    return py_minimize_molecules_back(mol1_pos, mol2_pos, mol2_com, sigma_scale,
                                      Nsteps, restrains, bonds_info, False,
                                      anchura, sim_types)

cdef std_matrix py_minimize_molecules_back(std_matrix mol1_pos,
                                           std_matrix mol2_pos,
                                           vector[double] mol2_com,
                                           double sigma_scale,
                                           long Nsteps,
                                           restrains_list restrains,
                                           bonds_list bonds_info,
                                           bool same_com,
                                           double anchura,
                                           vector[int] sim_types):

    return minimize_molecules(mol1_pos, mol2_pos, mol2_com, sigma_scale, Nsteps,
                              restrains, bonds_info, same_com, anchura,
                              sim_types)
