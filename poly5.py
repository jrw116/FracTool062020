# function to calculate the area of a polygon with 4 or more sides defined by nodes in 3D
# import functions used within this function
import numpy as np
import math

def poly5(point3,tol,g) :
    angles = {}
    # create reference vector
    ref = point3[1, :] - point3[0, :]
    angles[1] = 0
    # iterate through other points making up the polygon
    for l in range(2, len(point3)):
        ref2 = point3[l, :] - point3[0, :]
        cros = np.cross(ref, ref2)
        product = np.linalg.norm(ref) * np.linalg.norm(ref2)
        dot = np.dot(ref, ref2)
        cos = dot / product
        if abs(point3[0,g]-point3[1,g]) < tol and point3[l,g] > point3[0,g] :
            angles[l] = math.acos(cos)
        elif abs(point3[0,g]-point3[1,g]) < tol and point3[l,g] < point3[0,g] :
            angles[l] = -math.acos(cos)
        elif cros[2] < 0:
            angles[l] = -math.acos(cos)
        else:
            angles[l] = math.acos(cos)
    # create sorted angles array
    sorted_angles = sorted(angles.items(), key=lambda x: x[1])
    # create sorted index array of polygon points
    sorted_index = np.array(sorted_angles)[:, 0]
    area = 0
    # calculate area of the triangles making up the polygon
    for l in range(0, len(sorted_index) - 1):
        area2 = 0.5 * np.linalg.norm(
            np.cross((point3[sorted_index[l], :] - point3[0, :]), (point3[sorted_index[l + 1], :] - point3[0, :])))
        area = area + area2
    return(area)