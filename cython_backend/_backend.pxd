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
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.map cimport map as cppmap
from libcpp.string cimport string
from libcpp cimport bool

ctypedef vector[vector[double]] std_matrix
ctypedef vector[vector[int]] restrains_list
ctypedef cppmap[int, vector[pair[int, double]]] bonds_list

cdef extern from "c++_src/gaddle_backend.h":
    cdef std_matrix sqeuclidean_distance(std_matrix ma1,
                                         std_matrix ma2)
cdef extern from "c++_src/gaddle_backend.h":
    cdef std_matrix minimize_molecules(std_matrix mol1_pos,
                                       std_matrix mol2_pos,
                                       vector[double] mol2_com,
                                       double sigma_scale,
                                       long Nsteps,
                                       restrains_list restrains,
                                       bonds_list bonds_info,
                                       bool same_com,
                                       double anchura,
                                       vector[int]sim_types)
