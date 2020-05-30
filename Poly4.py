# function to calculate the area of a polygon with 4 sides defined by nodes in 3D
# import functions used within this function
import numpy as np
import math

def poly4(point3) :
    q = point3
    p1 = np.zeros((3, 4))
    # half each side of the quadrilateral
    p1[:, 0] = (q[0, :] + q[1, :]) / 2
    p1[:, 1] = (q[1, :] + q[2, :]) / 2
    p1[:, 2] = (q[2, :] + q[3, :]) / 2
    p1[:, 3] = (q[3, :] + q[0, :]) / 2
    cross = np.zeros((3, 1))
    cross[0] = (p1[1, 1] - p1[1, 0]) * (p1[2, 2] - p1[2, 0]) - (p1[2, 1] - p1[2, 0]) * (
            p1[1, 2] - p1[1, 0])
    cross[1] = (p1[2, 1] - p1[2, 0]) * (p1[0, 2] - p1[0, 0]) - (p1[0, 1] - p1[0, 0]) * (
            p1[2, 2] - p1[2, 0])
    cross[2] = (p1[0, 1] - p1[0, 0]) * (p1[1, 2] - p1[1, 0]) - (p1[1, 1] - p1[1, 0]) * (
            p1[0, 2] - p1[0, 0])
    # calculate area of diamond shape
    area = math.sqrt(np.sum(np.square(cross)))
    # double area to get area of quadrilateral
    area1 = area * 2
    # change sorting of the points
    point2 = point3
    point2[2, :], point2[3, :] = point2[3, :], point2[2, :].copy()

    q = point2
    p1 = np.zeros((3, 4))
    # half each side of the quadrilateral
    p1[:, 0] = (q[0, :] + q[1, :]) / 2
    p1[:, 1] = (q[1, :] + q[2, :]) / 2
    p1[:, 2] = (q[2, :] + q[3, :]) / 2
    p1[:, 3] = (q[3, :] + q[0, :]) / 2
    p1[:, 3] = (q[3, :] + q[0, :]) / 2
    cross = np.zeros((3, 1))
    cross[0] = (p1[1, 1] - p1[1, 0]) * (p1[2, 2] - p1[2, 0]) - (p1[2, 1] - p1[2, 0]) * (
            p1[1, 2] - p1[1, 0])
    cross[1] = (p1[2, 1] - p1[2, 0]) * (p1[0, 2] - p1[0, 0]) - (p1[0, 1] - p1[0, 0]) * (
            p1[2, 2] - p1[2, 0])
    cross[2] = (p1[0, 1] - p1[0, 0]) * (p1[1, 2] - p1[1, 0]) - (p1[1, 1] - p1[1, 0]) * (
            p1[0, 2] - p1[0, 0])
    # calculate area of diamond shape
    area = math.sqrt(np.sum(np.square(cross)))
    # double area to get area of quadrilateral
    area2 = area * 2
    # the sorting of the points that produces the greatest area is the correct sorting
    if area1 < area2:
        point3 = point2
        area = area2
    else:
        area = area1
    # output the area and the points sorted correctly
    return(area, point3)