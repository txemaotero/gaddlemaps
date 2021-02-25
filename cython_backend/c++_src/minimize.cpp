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
#include "minimize.h"
#include "algorithm_utils.h"
#include "transform_molecule.h"

const double PI = 3.141592653589793238462643383279502884;
const double PI4 = PI/4.;

arma::mat rotation_matrix(arma::mat &axis, const double theta){
    arma::mat eye(3,3);
    eye.eye();
    axis /= arma::norm(axis);

    arma::mat ddt = axis.t() * axis;
    arma::mat skew;
    skew << 0 << axis(2) << -axis(1) << arma::endr
              << -axis(2) << 0 << axis(0) << arma::endr
              << axis(1) << -axis(0) << 0 << arma::endr;


    return ddt + std::cos(theta) * (eye - ddt) + std::sin(theta) * skew;
}

auto random_engine3 = std::default_random_engine();
arma::mat minimize_molecules_routine(arma::mat &mol1_pos, arma::mat &mol2_pos,
                                arma:: rowvec mol2_com, double sigma_scale,
                                long Nsteps, restraint_list &restrains,
                                bonds &bonds_info, bool same_com, double anchura,
                                const std::vector<int> sim_types){
    long counter = 0;
    long global_counter = 0;
    auto rand_sim = std::uniform_int_distribution<int>(0,sim_types.size()-1);
    auto rand_theta = std::normal_distribution<double>(0, PI4);

    arma::mat test_pos = mol2_pos;
    arma::mat mol2_opt = mol2_pos;

    auto chi2_molecules = Chi2Calculator(mol1_pos, mol2_pos, restrains);

    double chi2 = chi2_molecules(mol2_pos);
    double chi2min = chi2;
    double chi2_new = 0;
    arma::rowvec desplazamiento(3);
    auto progressbar = Progressbar(Nsteps, chi2min);

    std::cout << "\n\n";

    while (counter < Nsteps){
        test_pos = mol2_pos;
        int change = sim_types[rand_sim(random_engine3)];

        // Rigid translations
        if (change == 0){
            if (same_com){
                continue;
            }
            desplazamiento = arma::randn<arma::rowvec>(3)*anchura;
            test_pos.each_row() += desplazamiento;
        }

        // Rigid rotations
        else if( change == 1){
            arma::rowvec axis = (arma::randu<arma::rowvec>(3)*2)-1;
            double theta = rand_theta(random_engine3);
            auto matrix_rot = rotation_matrix(axis, theta);
            test_pos = (mol2_pos.each_row()-mol2_com)*matrix_rot;
            test_pos.each_row() += mol2_com;
        }

        else if (change == 2){
            move_mol_atom(test_pos, bonds_info, sigma_scale);
        }

        // General changes common to all the movements
        chi2_new = chi2_molecules(test_pos);
        if (accept_metropolis(chi2, chi2_new)){
            mol2_pos = test_pos;
            chi2 = chi2_new;

            // Correct the position of center of mass if needed
            if (change == 0){
                mol2_com += desplazamiento;
            }
        }

        if (chi2 < chi2min){
            chi2min = chi2;
            mol2_opt = mol2_pos;
            // std::cout << "\r'\tChi2 =" << chi2min
            //           << " " << global_counter  << " " << counter << std::flush;
            counter = -1;

        }
        counter++;
        global_counter ++;
        progressbar(counter, chi2min);
    }
    std::cout << "\n" << counter << " " << global_counter << std::endl;
    return mol2_opt;
}
