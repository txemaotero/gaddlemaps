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
"""
This submodule contains some useful functions that are used in other parts of
the package.
"""

from typing import List, Tuple

import numpy as np


def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    """
    Returns the rotation matrix associated to an angle and an axis.

    Parameters
    ----------
    axis : list or numpy.ndarray
        A vector in 3D space that defines the axis of rotation
    theta : float
        The angle to rotate in radians in counter clockwise.

    """
    norm_ax = axis / np.linalg.norm(axis)

    eye = np.eye(3, dtype=np.float64)
    ddt = np.outer(norm_ax, norm_ax)
    skew = np.array([[0, norm_ax[2], -norm_ax[1]],
                     [-norm_ax[2], 0, norm_ax[0]],
                     [norm_ax[1], -norm_ax[0], 0]], dtype=np.float64)
    mtx = ddt + np.cos(theta) * (eye - ddt) + np.sin(theta) * skew
    return mtx


def calcule_base(pos: List[np.ndarray]) -> Tuple[Tuple[np.ndarray, ...],
                                                 np.ndarray]:
    """
    Calculates a orthonormal base from a list of 3 atoms.

    Given a list of three vectors with the position of three atoms, this
    function returns a vector basis and the application point. The first
    vector goes from the application point to the last atom. The second one
    is normal to the plane which contains the three atoms and perpendicular
    to the first vector. The last is perpendicular to the others. All are
    unitary forming an orthonormal basis. In case of collinear points, the
    second vector is set to ([vec1[1], -vec1[0], 0]).

    Parameters
    ----------
    pos : list
        List with three vectors (numpy.ndarray) with the position of the
        atoms.

    Returns
    -------
    base : tuple of numpy.ndarray
        A tuple with the three vectors of the base.
    app_point : float
        The application point of the base. It corresponds to pos[0]

    """

    pos0, pos1, pos2 = pos
    vec1 = pos2-pos0
    vec1 /= np.linalg.norm(vec1)
    vec3 = np.cross(vec1, pos1-pos0)
    # Consider the aligned vector case
    if not np.any(vec3):
        v10, v11, _ = vec1
        vec3 = np.array([v11, -v10, 0]) / (v10**2+v11**2)
    else:
        vec3 /= np.linalg.norm(vec3)
    vec2 = np.cross(vec3, vec1)
    return (vec1, vec2, vec3), pos0
