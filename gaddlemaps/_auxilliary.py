import numpy as np


def rotation_matrix(axis, theta):
    """
    Returns the rotation matrix associated to an angle with cos, sin
    values of the cosine and sine.

    Parameters
    ----------
    axis : list or numpy.ndarray
        A vector that defines the axis of rotation
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


def calcule_base(pos):
    """
    Calculates a orthonormal base from a list of 3 atoms.

    Given a list of three vectors with the position of three atoms, this
    function returns a vector basis and the application point. The first
    vector goes from the application point to the last specified atom. The
    second one is normal to the plane which contains the three atoms and
    perpendicular to the first vector. The last is perpendicular to the
    others. All are unitary forming an ortonormal basis. In case of colinear
    points, the second vector is set to ([vec1[1], -vec1[0], 0]).

    Parameters:
    -----------
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
