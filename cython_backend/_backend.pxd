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
