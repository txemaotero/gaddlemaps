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
