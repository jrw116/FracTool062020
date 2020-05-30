# function to calculate the cut plane area within each tetrahedral element of a meshed model
# import functions used within this function
from Poly4 import poly4
from Poly5 import poly5
from Line import line
import numpy as np
import math

def plycomptetel(e, nodeSet, p, perp, t, tol, area3, areaelply, Vert, minVert, r, f, g, h, x, perp1, nPlies, m, sympos, symneg, m1) :
    # pull node connectivity from model information
    connected = p.elements[e].connectivity
    sort = np.sort(connected)
    label = p.elements[e].label
    # workout length of sides in element
    side1 = np.array(nodeSet[sort[0]].coordinates) - np.array(nodeSet[sort[1]].coordinates)
    side2 = np.array(nodeSet[sort[0]].coordinates) - np.array(nodeSet[sort[2]].coordinates)
    side3 = np.array(nodeSet[sort[0]].coordinates) - np.array(nodeSet[sort[3]].coordinates)
    size = [np.linalg.norm(side1)]
    size.append(np.linalg.norm(side2))
    size.append(np.linalg.norm(side3))
    maxSize = max(size)
    diag = math.sqrt(maxSize ** 2 + maxSize ** 2)
    # check if elements are close to plane and if they are calculate the area
    if float(abs(np.dot(perp, np.array(nodeSet[sort[0]].coordinates)) - m1) / (np.linalg.norm(perp))) <= float(
            math.sqrt(diag ** 2 + maxSize ** 2)):
        # initialise variables
        vertex = np.zeros((4, 3))
        point = np.zeros((4, 3))
        area3[label] = []
        areaelply[label] = [0] * nPlies
        # assign vertices of element
        vertex[0, :] = np.array(nodeSet[sort[0]].coordinates)
        vertex[1, :] = np.array(nodeSet[sort[1]].coordinates)
        vertex[2, :] = np.array(nodeSet[sort[2]].coordinates)
        vertex[3, :] = np.array(nodeSet[sort[3]].coordinates)
        # find position of each node of the element with respect to the cut plane
        d = np.zeros((4, 1))
        for i in range(0, 4):
            d[i] = -np.dot(perp, vertex[i, :]) + np.dot(perp, t)
        # set variable value to latter check the cut plane intersects the model
        for i in range(0, len(d)):
            if d[i] > 0:
                sympos = 1
            elif d[i] < 0:
                symneg = 1
        # iterate through each node to see if any lie on the cut plane
        l = 0
        for j in range(0, 4):
            if d[j] < 0 or d[j] > 0:
                l = l + 1
            else:
                l = l
        if l != 0:
            intersection = 0
        j3 = 0
        for j1 in range(0, 4):
            # if node lies on cut plane assign to point
            if d[j1] == 0:
                intersection += 1
                point[intersection - 1, :] = vertex[j1, :]
                j3 += 1
            # if two nodes lie either side of the cut plane find the intersection point of the edge
            # between the two nodes and assign to point
            else:
                for j2 in range(j1 + 1, 4):
                    if d[j1] * d[j2] < 0:
                        intersection += 1
                        point[intersection - 1, :] = (d[j1] * vertex[j2, :] - d[j2] * vertex[j1, :]) / (
                                d[j1] - d[j2])
        # if there are three intersection points between the cut plane and the element find the area
        if intersection == 3:
            point = np.delete(point, 3, 0)
            area = 0.5 * np.linalg.norm(np.cross((point[1, :] - point[0, :]), (point[2, :] - point[0, :])))
        # if there are four intersection points between the cut plane and the element find the area
        elif intersection == 4:
            area, point = poly4(point)
            m += 1
        # if there are no intersection points set the area to zero
        else:
            area = 0
        area3[label].append(area)

        # if there is cut plane area within the element enter this loop to split the area up by the amount in each ply
        if intersection == 3 or intersection == 4:
            # set counters
            t1 = 0
            t3 = 0
            t4 = 0
            t5 = 0
            # iterate through each ply
            for z in range(0, len(Vert)):
                # if the intersection points inside the ply assign previously calculated area
                if np.all(point[:, r] >= Vert[z - 1] - tol) and np.all(point[:, r] <= Vert[z] + tol) and z == (
                        len(Vert) - 1) and t5 == 0:
                    areaelply[label][z] = area3[label][0]
                    t5 = 1
                # if the intersection points lie outside the ply set the area to zero for that ply
                elif np.all(point[:, r] >= Vert[z - 1] - tol) and np.all(
                        point[:, r] <= Vert[z] + tol) and z != 0 and t5 == 0:
                    areaelply[label][z] = area3[label][0]
                    t5 = 1
                # if the intersection points lie outside the ply set the area to zero for that ply
                elif np.all(point[:, r] >= minVert - tol) and np.all(
                        point[:, r] <= Vert[z] + tol) and z == 0 and t5 == 0:
                    areaelply[label][z] = area3[label][0]
                    t5 = 1
                # if all the intersection points lie outside the ply skip this iteration
                elif np.all(point[:, r] > Vert[z] - tol):
                    continue
                # if all the intersection points lie outside the ply skip this iteration
                elif np.all(point[:, r] < Vert[z - 1] + tol) and z != 0:
                    continue
                # if the intersection points inside the ply enter this loop
                elif t4 != 1:
                    # initialise variables
                    point1 = np.zeros((2, 3))
                    d1 = np.zeros((intersection, 1))
                    # find position of each node of the element with respect to the ply top surface/plane
                    for i in range(0, intersection):
                        d1[i] = -np.dot(perp1[z, :], point[i, :]) + np.dot(perp1[z, :], x[z, :])
                    # if there are three intersection points enter this loop to find the points of intersection
                    # between the ply top surface and the cut plane area within the element
                    if intersection == 3:
                        l = 0
                        for j in range(0, intersection):
                            if d1[j] < 0 or d1[j] > 0:
                                l = l + 1
                            else:
                                l = l
                        if l != 0:
                            intersection1 = 0
                        for j1 in range(0, intersection):
                            if d1[j1] == 0:
                                intersection1 += 1
                                point1[intersection1 - 1, :] = point[j1, :]
                            else:
                                for j2 in range(j1 + 1, intersection):
                                    if d1[j1] * d1[j2] < 0:
                                        intersection1 += 1
                                        point1[intersection1 - 1, :] = (d1[j1] * point[j2, :] - d1[j2] * point[j1,
                                                                                                         :]) / (
                                                                               d1[j1] - d1[j2])
                    # if there are four intersection points enter this loop to find the points of intersection
                    # between the ply top surface and the cut plane area within the element
                    if intersection == 4:
                        l = 0
                        for j in range(0, intersection):
                            if d1[j] < 0 or d1[j] > 0:
                                l = l + 1
                            else:
                                l = l
                        if l != 0:
                            intersection1 = 0
                        for j1 in range(0, intersection - 1):
                            if d1[j1] == 0:
                                intersection1 += 1
                                point1[intersection1 - 1, :] = point[j1, :]
                            else:
                                if d1[j1] * d1[j1 + 1] < 0:
                                    intersection1 += 1
                                    point1[intersection1 - 1, :] = (d1[j1] * point[j1 + 1, :] - d1[j1 + 1] * point[j1,
                                                                                                             :]) / (
                                                                           d1[j1] - d1[j1 + 1])
                        # loop round to compare 4th and 1st intersection points between the cut plane and element edges
                        if d1[intersection - 1] == 0:
                            intersection1 += 1
                            point1[intersection1 - 1, :] = point[intersection - 1, :]
                        else:
                            if d1[intersection - 1] * d1[0] < 0:
                                intersection1 += 1
                                point1[intersection1 - 1, :] = (d1[intersection - 1] * point[0, :] - d1[0] * point[
                                                                                                             intersection - 1,
                                                                                                             :]) / (
                                                                       d1[intersection - 1] - d1[0])
                    # if no area has been found in a ply yet enter this loop
                    if t1 == 0:
                        # initialise counter and variable
                        k = 0
                        pointbelow = np.zeros((1, 3))
                        # run through each intersection point between the element edges and the cut plane
                        for i in range(0, intersection):
                            # if intersection point lies between the top and bottom surfaces of the ply assign to pointbelow
                            if point[i, r] < Vert[z]:
                                if k == 0:
                                    pointbelow[k, :] = point[i, :].copy()
                                if k != 0:
                                    pointbelow = np.vstack((pointbelow, point[i, :]))
                                k += 1

                        # combine pointbelow and intersection points between the ply top surface and the eges of the cut plane within the element
                        point3 = np.vstack((point1, pointbelow))
                        t1 += 1
                        # if the shape of the area has 3 nodes find area of the triangle
                        if len(point3) == 3:
                            area6 = 0.5 * np.linalg.norm(
                                np.cross((point3[1, :] - point3[0, :]), (point3[2, :] - point3[0, :])))
                            areaelply[label][z] = areaelply[label][z] + area6
                        # if the shape of the area has 4 nodes find area of the square using poly4 function
                        if len(point3) == 4:
                            area, point3 = poly4(point3)
                            areaelply[label][z] = areaelply[label][z] + area
                        # if the shape of the area has 5 nodes find area of the pentagon using poly5 function
                        if len(point3) == 5:
                            area = poly5(point3, tol, g)
                            areaelply[label][z] = areaelply[label][z] + area
                    # enter this loop if the first intersection point between plies has been found
                    elif z != len(Vert) - 1 and t1 != 0:
                        # initialise variables and counter
                        point3 = np.vstack((point4, point1))
                        pointother = np.empty([1, 3])
                        t1 += 1
                        # iterate through all the intersection points
                        for i in range(0, intersection):
                            # initialise counters
                            count = 0
                            count1 = 0
                            count2 = 0
                            count3 = 0
                            count4 = 0
                            # compare each intersection point to all other intersection points
                            for j in range(0, intersection):
                                # only enter loop if no jth intersection point has been found that lies in the same ply as the ith intersection point
                                if j != i and np.all(point[i, :] != pointother):
                                    # if both the ith and jth intersection point lie in the plane assign the jth intersection point to pointother
                                    if Vert[z - 1] < point[j, r] < Vert[z] and Vert[z - 1] < point[i, r] < Vert[z]:
                                        count4 = 1
                                        pointother = point[j, :].copy()
                            # for the 1st and 2nd points if the number of intersections = 3 or for the 1st, 2nd and 3rd points if the number of intersections = 4
                            if i <= intersection - 2:
                                # if the point lies in the ply and there are no other points lying in the ply enter this loop
                                if Vert[z - 1] < point[i, r] < Vert[z] and t3 != 1 and count4 != 1:
                                    point5 = np.zeros((3, 3))
                                    # run the function line to determine which intersection points between the ply bottom and top surfaces
                                    # and the cut plane area shape in the element are closest to the node of the cut plane area in question
                                    point5 = line(point5, point3, point, tol, i, f, g, h)
                                    point5[2, :] = point[i, :]
                                    # calculate the area of the triangle
                                    area9 = 0.5 * np.linalg.norm(
                                        np.cross((point5[1, :] - point5[0, :]), (point5[2, :] - point5[0, :])))
                                    areaelply[label][z] = areaelply[label][z] + area9
                                # if there are 4 intersection points between the element and the cut plane enter this loop
                                elif intersection == 4:
                                    # if the intersection point lies in the ply and it is known there is another intersection point also lying in the ply then enter this loop
                                    if Vert[z - 1] < point[i, r] < Vert[z] and count4 == 1:
                                        for l in range(0, intersection):
                                            # depending on the stack direction if the x, y or z values of the first 3 points are the same enter this loop
                                            if abs(point[0, f] - point[1, f]) < tol and abs(
                                                    point[1, f] - point[2, f]) < tol:
                                                # depending on stack direction create counters depending on intersection point position compared to position of intersection
                                                # points between ply top and bottom surfaces and the edges of the cut plane area within the element
                                                if point[i, h] > point3[l, h]:
                                                    count = count + 1
                                                if point[i, h] < point3[l, h]:
                                                    count1 = count1 + 1
                                                if pointother[h] > point3[l, h]:
                                                    count2 = count2 + 1
                                                if pointother[h] < point3[l, h]:
                                                    count3 = count3 + 1
                                            else:
                                                if point[i, f] > point3[l, f]:
                                                    count = count + 1
                                                if point[i, f] < point3[l, f]:
                                                    count1 = count1 + 1
                                                if pointother[f] > point3[l, f]:
                                                    count2 = count2 + 1
                                                if pointother[f] < point3[l, f]:
                                                    count3 = count3 + 1
                                        # if both points are both to the left or both to the right of the ply surface intersection points enter this loop
                                        if (count >= 3 and count2 >= 3) or (count <= 1 and count2 <= 1):
                                            point5 = np.zeros((4, 3))
                                            # run the function line to determine which intersection points between the ply bottom and top surfaces
                                            # and the cut plane area shape in the element are closest to the node of the cut plane area in question
                                            point5 = line(point5, point3, point, tol, i, f, g, h)
                                            point5[2, :] = point[i, :].copy()
                                            point5[3, :] = pointother.copy()
                                            # calculate the area of the quadrilateral
                                            area, point5 = poly4(point5)
                                            areaelply[label][z] = areaelply[label][z] + area
                                            t3 = 1
                                        # if only one intersection point in the ply enter this loop
                                        else:
                                            point5 = np.zeros((3, 3))
                                            # run the function line to determine which intersection points between the ply bottom and top surfaces
                                            # and the cut plane area shape in the element are closest to the node of the cut plane area in question
                                            point5 = line(point5, point3, point, tol, i, f, g, h)
                                            point5[2, :] = point[i, :]
                                            # calculate the area of the triangle
                                            area9 = 0.5 * np.linalg.norm(
                                                np.cross((point5[1, :] - point5[0, :]), (point5[2, :] - point5[0, :])))
                                            areaelply[label][z] = areaelply[label][z] + area9
                            # for the 3rd point if the number of intersections = 3 or for the 4th point if the number of intersections = 4
                            elif i == intersection - 1:
                                # if the intersection point is the only one lying in the ply enter this loop
                                if Vert[z - 1] < point[i, r] < Vert[z] and t3 != 1 and count4 != 1:
                                    point5 = np.zeros((3, 3))
                                    # run the function line to determine which intersection points between the ply bottom and top surfaces
                                    # and the cut plane area shape in the element are closest to the node of the cut plane area in question
                                    point5 = line(point5, point3, point, tol, i, f, g, h)
                                    point5[2, :] = point[i, :]
                                    # calculate the area of the triangle
                                    area9 = 0.5 * np.linalg.norm(
                                        np.cross((point5[1, :] - point5[0, :]), (point5[2, :] - point5[0, :])))
                                    areaelply[label][z] = areaelply[label][z] + area9
                        # calculate central area between all the intersection points of the top and bottom surfaces of the ply and the cut plane area within the element
                        if len(point3) == 4:
                            area, point3 = poly4(point3)
                            areaelply[label][z] = areaelply[label][z] + area
                    # enter this loop as long as not on the last ply
                    if z != len(Vert) - 1:
                        # enter this loop if all intersection points are less than the top surface of the final ply
                        if np.all(point[:, r] <= (Vert[z + 1] + tol)):
                            # initialise counters and variables
                            t4 = 1
                            k = 0
                            pointabove = np.zeros((1, 3))
                            # iterate through intersection points
                            for i in range(0, intersection):
                                # if intersection point is in the final ply enter this loop
                                if point[i, r] > Vert[z]:
                                    if k == 0:
                                        pointabove[k, :] = point[i, :].copy()
                                    if k != 0:
                                        pointabove = np.vstack((pointabove, point[i, :]))
                                    k += 1
                            # combine pointabove and intersection points between the ply top surface and the eges of the cut plane within the element
                            point3 = np.vstack((point1, pointabove))
                            # if the area shape produced is a triangle find the area
                            if len(point3) == 3:
                                area6 = 0.5 * np.linalg.norm(
                                    np.cross((point3[1, :] - point3[0, :]), (point3[2, :] - point3[0, :])))
                                areaelply[label][z + 1] = areaelply[label][z + 1] + area6
                            # if the area shape produced is a quadrilateral find the area
                            if len(point3) == 4:
                                area, point3 = poly4(point3)
                                areaelply[label][z + 1] = areaelply[label][z + 1] + area
                            # if the area shape produced is a pentagon find the area
                            if len(point3) == 5:
                                area = poly5(point3, tol, g)
                                areaelply[label][z + 1] = areaelply[label][z + 1] + area
                    # enter this loop if on the final ply
                    if z == len(Vert) - 1:
                        # enter this loop if any of the intersection points are on the last ply's top surface
                        if np.any(abs(point[:, r] - Vert[z]) < tol):
                            # initialise counters and variables
                            t4 = 1
                            k = 0
                            pointabove = np.zeros((1, 3))
                            # iterate through all points to see if any of the intersection points lie in the final ply or on the top surface of the final ply
                            for i in range(0, intersection):
                                if abs(point[i, r] - Vert[z]) < tol or Vert[z - 1] < point[i, r] < Vert[z]:
                                    if k == 0:
                                        pointabove[k, :] = point[i, :].copy()
                                    if k != 0:
                                        pointabove = np.vstack((pointabove, point[i, :]))
                                    k += 1
                            # combine pointabove and intersection points between the ply top surface and the eges of
                            # the cut plane within the element
                            point3 = np.vstack((point4, pointabove))
                            # if the area shape produced is a triangle find the area
                            if len(point3) == 3:
                                area6 = 0.5 * np.linalg.norm(
                                    np.cross((point3[1, :] - point3[0, :]), (point3[2, :] - point3[0, :])))
                                areaelply[label][z] = areaelply[label][z] + area6
                            # if the area shape produced is a quadrilateral find the area
                            if len(point3) == 4:
                                area, point3 = poly4(point3)
                                areaelply[label][z] = areaelply[label][z] + area
                            # if the area shape produced is a pentagon find the area
                            if len(point3) == 5:
                                area = poly5(point3, tol, g)
                                areaelply[label][z] = areaelply[label][z] + area
                    # set top intersection points between the ply top surface and the cut plane area in the shape to
                    # point4 for use with next loop for next ply
                    point4 = point1
        # if the cut plane lies on the face of the element half the area on that cut face
        if j3 == 3:
            for j in range(0, nPlies):
                areaelply[label][j] = 0.5 * areaelply[label][j]
        areaelply[label].reverse()
    # if element does not lie close to cut plane then the cut plane area is 0
    else :
        areaelply[label] = [0] * nPlies
    return (areaelply, sympos, symneg)