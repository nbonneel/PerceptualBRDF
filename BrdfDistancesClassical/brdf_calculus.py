import math
import numpy as np
#import scipy.linalg

def get_intensity(colorRGB):
    # ITU-R BT.709
    #y_coef = [0.2125, 0.7154, 0.0721]    
    y_coef = [0.299, 0.587, 0.114]
    return sum(color*coef for color, coef in zip(colorRGB, y_coef))

def cartesian_to_polar(point):
    return math.sqrt(point[0] * point[0] + point[1] * point[1]) , math.atan2(point[1], point[0])

def polar_to_cartesian(point):
    return point[0] * math.cos(point[1]), point[0] * math.sin(point[1])

def spherical_to_cartesian(theta, phi):
    proj_vec = np.sin(theta)
    
    vec = np.array([proj_vec * np.cos(phi),
                    proj_vec * np.sin(phi),
                    np.cos(theta)])

    return vec # / np.linalg.norm(vec)

def cartesian_to_spherical(vector):
    """
    :return: (theta, phi)
    """
    v = vector / np.linalg.norm(vector)
    return np.arccos(v[2]), np.arctan2(v[1], v[0])
    
# https://fr.wikipedia.org/wiki/Aire_et_centre_de_masse_d%27un_polygone
def get_gravity_center_polar(points_polar):
    gravity = cartesian_to_polar(get_gravity_center_cartesian(points_polar))

    if gravity[0] < 0:
        gravity =  -gravity[0], gravity[1] + math.pi
        
    return gravity

def get_gravity_center_cartesian(points_polar):
    points_cartesian = [polar_to_cartesian(p) for p in points_polar]
    area = get_area(points_cartesian)

    gravity = [0]*2

    for i in range(0, len(points_cartesian) - 1):
        coef = points_cartesian[i][0] * points_cartesian[i+1][1] - points_cartesian[i+1][0] * points_cartesian[i][1]
        gravity[0] += (points_cartesian[i][0] + points_cartesian[i+1][0]) * coef
        gravity[1] += (points_cartesian[i][1] + points_cartesian[i+1][1]) * coef
        
    gravity[0] /= (6 * area)
    gravity[1] /= (6 * area)

    return gravity

def get_area(points_cartesian):
    a = 0
    for i in range(0, len(points_cartesian) - 1):
        a += points_cartesian[i][0] * points_cartesian[i+1][1] - points_cartesian[i+1][0] * points_cartesian[i][1]

    return a/2

#def rotation_matrix2(axis, theta):
#    return scipy.linalg.expm3(np.cross(np.eye(3), axis/np.linalg.norm(axis)*theta))

# http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

def rotate_vector(vector, axis, angle):
    cos_ang = math.cos(angle)
    sin_ang = math.sin(angle)

    out = vector * cos_ang
    tmp = np.dot(axis, vector) * (1 - cos_ang)

    out += np.array(axis) * tmp
    cross = np.cross(axis, vector)
    out += cross * sin_ang
    return out

def std_coords_to_half_diff(theta_i, phi_i, theta_o, phi_o):
    in_vect = spherical_to_cartesian(theta_i, phi_i)
    out_vect = spherical_to_cartesian(theta_o, phi_o)
    
    half_vect = in_vect + out_vect

    theta_h, phi_h = cartesian_to_spherical(half_vect)

    diff_vect = np.matmul(rotation_matrix([0, 0, 1], -phi_h), in_vect)
    diff_vect = np.matmul(rotation_matrix([0, 1, 0], -theta_h), diff_vect)

    theta_d, phi_d = cartesian_to_spherical(diff_vect)

    return theta_h, phi_h, theta_d, phi_d

def half_diff_coords_to_std(theta_h, phi_h, theta_d, phi_d):
    diff_vect = spherical_to_cartesian(theta_d, phi_d)
    half_vect = spherical_to_cartesian(theta_h, phi_h)

    in_vect = np.matmul(rotation_matrix([0, 1, 0], theta_h), diff_vect)
    in_vect = np.matmul(rotation_matrix([0, 0, 1], phi_h), in_vect)

    in_vect /= np.linalg.norm(in_vect)

    theta_i, phi_i = cartesian_to_spherical(in_vect)

    l = np.dot(half_vect, in_vect)
    out_vect = 2 * l * half_vect - in_vect
    
    theta_o, phi_o = cartesian_to_spherical(out_vect)
    
    return theta_i, phi_i, theta_o, phi_o
    
def get_vector_after_ns(vector, theta_n, phi_n):
    # Verbose:    
    # t_t_mat = calc.rotation_matrix(t, theta)
    # b_t_mat = calc.rotation_matrix(b, theta)
    # n_t_mat = calc.rotation_matrix(n, theta)
    
    # t_p_mat = calc.rotation_matrix(t, phi)
    # b_p_mat = calc.rotation_matrix(b, phi)
    # n_p_mat = calc.rotation_matrix(n, phi)

    # t_s2 = np.matmul(n_t_mat, t)
    # t_s2 = np.matmul(t_p_mat, t_s2)
    
    # b_s2 = np.matmul(t_t_mat, b)
    # b_s2 = np.matmul(b_p_mat, b_s2)
    
    # n_s2 = np.matmul(b_t_mat, n)
    # n_s2 = np.matmul(n_p_mat, n_s2)

    n_s = spherical_to_cartesian(theta_n, phi_n)
    t_s = np.array([n_s[2], n_s[0], n_s[1]])
    b_s = np.array([n_s[1], n_s[2], n_s[0]])
    tbn = np.transpose([t_s, b_s, n_s])

    return np.matmul(tbn, vector)

def get_vector_spherical_after_ns(theta, phi, theta_n, phi_n):
    vec = spherical_to_cartesian(theta, phi)
    vec = get_vector_after_ns(vec, theta_n, phi_n)
    
    return cartesian_to_spherical(vec)
