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
#include "transform_molecule.h"
auto random_engine2 = std::default_random_engine();
auto random_direction = std::uniform_int_distribution<int>(0,1.);
auto random_normal = std::normal_distribution<double>(0.,1.);

void move_mol_atom(arma::mat &atom_pos,
                   bonds &bonds_info,double sigma_scale){
    int Natoms = int(atom_pos.n_rows);
    auto int_dist = std::uniform_int_distribution<int>(0, Natoms-1);
    // Choose one random atom to move
    int atom_index = int_dist(random_engine2);
    // Choose a random displacement
    auto disp = find_atom_random_displ(atom_pos, bonds_info,
                                       atom_index, sigma_scale);

    // Start an empty set of the atoms to move
    std::set<int> to_move;
    // Poulate the queue
    for (int atom=0; atom < Natoms; ++atom){
        to_move.insert(atom);
    }
    // Move the first atom
    atom_pos.row(atom_index) += disp;
    // Remove the value from the queue
    to_move.erase(atom_index);

    // Create the waiting list of the atoms to move
    std::queue<std::pair<int, std::pair<int, double>>> waiting;

    auto bonds_list = bonds_info[atom_index];
    for (int i=0; i < bonds_list.size(); ++i){
        waiting.push(std::make_pair(atom_index, bonds_list[i]));
        to_move.erase((bonds_list[atom_index].first));
    }

    while (! waiting.empty()){
        // Get the info of the priority atom to move
        auto bond_move = waiting.front();
        int ind1 = bond_move.first;
        int ind2 = bond_move.second.first;
        double bond = bond_move.second.second;
        waiting.pop();

        // Move the atom to fix the distance
        arma::rowvec diferencia = atom_pos.row(ind1) - atom_pos.row(ind2);
        double modulo = arma::norm(diferencia);
        diferencia /= modulo;
        atom_pos.row(ind2) += (modulo - bond)*diferencia;

        // Adds the atom bonded to the extreme atom to the queue and remove them
        // from the to_move list.
        bonds_list = bonds_info[ind2];
        for (int i=0; i < bonds_list.size(); ++i){
            auto new_bond = bonds_list[i].first;
            if (to_move.find(new_bond) != to_move.end()){
                waiting.push(std::make_pair(ind2, bonds_list[i]));
                to_move.erase(new_bond);
            }
        }
    }
}


arma::rowvec find_atom_random_displ(const arma::mat &atom_pos,
                                    bonds &bonds_info, int atom_index,
                                    double sigma_scale){
    int Nbonded_ref = int(bonds_info[atom_index].size());

    double sigma = bonds_info[atom_index][0].second * sigma_scale;

    arma::rowvec direction(3);
    arma::rowvec dir_vec1(3);
    arma::rowvec dir_vec2(3);
    switch (Nbonded_ref) {
    case 1:
        dir_vec1 = arma::randu<arma::rowvec>(3);
        dir_vec2 = atom_pos.row(bonds_info[atom_index][0].first);
        dir_vec2 -= atom_pos.row(atom_index);
        break;
    case 2:
        dir_vec1 = arma::randu<arma::rowvec>(3);
        dir_vec2 = atom_pos.row(bonds_info[atom_index][0].first);
        dir_vec2 -= atom_pos.row(bonds_info[atom_index][1].first);
        break;
    default:
        dir_vec1 = atom_pos.row(bonds_info[atom_index][0].first);
        dir_vec1 -= atom_pos.row(bonds_info[atom_index][2].first);
        dir_vec2 = atom_pos.row(bonds_info[atom_index][0].first);
        dir_vec2 -= atom_pos.row(bonds_info[atom_index][1].first);
        break;
    }
    direction = arma::cross(dir_vec1, dir_vec2);
    direction /= arma::norm(direction);
    direction *= (random_direction(random_engine2)*2-1);
    direction *= sigma*random_normal(random_engine2);
    return direction;
}
