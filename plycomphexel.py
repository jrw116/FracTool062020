# function to calculate the cut plane area within each hexahedral element of a meshed model
# import functions used within this function
from collections import Counter
from poly4 import poly4
from poly5 import poly5
import numpy as np
import math

def plycomphexel(e, nodeSet, perp, area2, plyElements, stackAxis, nPlies, m1, minVert, Vert, tol, vertexsubply, perp1, x, dir, t, sympos, symneg):
    r = t
    connected = plyElements[e].connectivity
    label = plyElements[e].label
    area2[label] = []
    # sort node index to keep track of lines and vertex
    sort = np.sort(connected)
    # workout length of sides in element
    side1 = np.array(nodeSet[sort[0]].coordinates) - np.array(nodeSet[sort[1]].coordinates)
    side2 = np.array(nodeSet[sort[0]].coordinates) - np.array(nodeSet[sort[2]].coordinates)
    side3 = np.array(nodeSet[sort[0]].coordinates) - np.array(nodeSet[sort[4]].coordinates)
    size = [np.linalg.norm(side1)]
    size.append(np.linalg.norm(side2))
    size.append(np.linalg.norm(side3))
    maxSize = max(size)
    diag = math.sqrt(maxSize ** 2 + maxSize ** 2)
    # check if elements are close to plane
    if float(abs(np.dot(perp, np.array(nodeSet[sort[0]].coordinates)) - m1) / (np.linalg.norm(perp))) <= float(
            math.sqrt(diag ** 2 + maxSize ** 2)):
        # Three possible directions for the ply stack
        # Vertex change depending on stack direction
        vertexhigh = np.zeros((1, 3))
        vertexlow = np.zeros((1, 3))
        point = np.zeros((12, 3))
        if str(stackAxis) == 'STACK_1':
            # sort element nodes by z coordinate
            vertextemp = np.zeros((8, 3))
            for i in range(0, 8):
                vertextemp[i, :] = np.array(nodeSet[sort[i]].coordinates)
            vertextemp = vertextemp[vertextemp[:, 2].argsort()]

            # split element nodes by z coordinate then sort by x coordinate then split again
            vertexint = vertextemp[0:4]
            vertexint = vertexint[vertexint[:, 0].argsort()]
            vertexlow = np.vstack((vertexlow, vertexint[0:2]))
            vertexhigh = np.vstack((vertexhigh, vertexint[2:4]))
            vertexint = vertextemp[4:8]
            vertexint = vertexint[vertexint[:, 0].argsort()]
            vertexlow = np.vstack((vertexlow, vertexint[0:2]))
            vertexhigh = np.vstack((vertexhigh, vertexint[2:4]))

            vertexlow = np.delete(vertexlow, 0, 0)
            vertexhigh = np.delete(vertexhigh, 0, 0)

            # sort in z again
            vertexlow = vertexlow[vertexlow[:, 2].argsort()]
            vertexhigh = vertexhigh[vertexhigh[:, 2].argsort()]

            # create rotational sorting in y
            vertexlow1 = vertexlow[0:2, :]
            vertexlow2 = vertexlow[2:4, :]
            vertexlow2 = vertexlow2[vertexlow2[:, 1].argsort()]
            vertexlow1 = vertexlow1[vertexlow1[:, 1].argsort()[::-1]]
            vertexlow = np.vstack((vertexlow1, vertexlow2))
            vertexhigh1 = vertexhigh[0:2, :]
            vertexhigh2 = vertexhigh[2:4, :]
            vertexhigh2 = vertexhigh2[vertexhigh2[:, 1].argsort()]
            vertexhigh1 = vertexhigh1[vertexhigh1[:, 1].argsort()[::-1]]
            vertexhigh = np.vstack((vertexhigh1, vertexhigh2))

            vertexsubply = np.vstack((vertexlow, vertexhigh))
        elif str(stackAxis) == 'STACK_2':
            # sort element nodes by x coordinate
            vertextemp = np.zeros((8, 3))
            for i in range(0, 8):
                vertextemp[i, :] = np.array(nodeSet[sort[i]].coordinates)
            vertextemp = vertextemp[vertextemp[:, 0].argsort()]

            # split element nodes by x coordinate then sort by y coordinate then split again
            vertexint = vertextemp[0:4]
            vertexint = vertexint[vertexint[:, 1].argsort()]
            vertexlow = np.vstack((vertexlow, vertexint[0:2]))
            vertexhigh = np.vstack((vertexhigh, vertexint[2:4]))
            vertexint = vertextemp[4:8]
            vertexint = vertexint[vertexint[:, 1].argsort()]
            vertexlow = np.vstack((vertexlow, vertexint[0:2]))
            vertexhigh = np.vstack((vertexhigh, vertexint[2:4]))

            vertexlow = np.delete(vertexlow, 0, 0)
            vertexhigh = np.delete(vertexhigh, 0, 0)

            # sort in x again
            vertexlow = vertexlow[vertexlow[:, 0].argsort()]
            vertexhigh = vertexhigh[vertexhigh[:, 0].argsort()]

            # create rotational sorting in z
            vertexlow1 = vertexlow[0:2, :]
            vertexlow2 = vertexlow[2:4, :]
            vertexlow2 = vertexlow2[vertexlow2[:, 2].argsort()]
            vertexlow1 = vertexlow1[vertexlow1[:, 2].argsort()[::-1]]
            vertexlow = np.vstack((vertexlow1, vertexlow2))
            vertexhigh1 = vertexhigh[0:2, :]
            vertexhigh2 = vertexhigh[2:4, :]
            vertexhigh2 = vertexhigh2[vertexhigh2[:, 2].argsort()]
            vertexhigh1 = vertexhigh1[vertexhigh1[:, 2].argsort()[::-1]]
            vertexhigh = np.vstack((vertexhigh1, vertexhigh2))

            vertexsubply = np.vstack((vertexlow, vertexhigh))
        elif str(stackAxis) == 'STACK_3':
            # sort element nodes by y coordinate
            vertextemp = np.zeros((8, 3))
            for i in range(0, 8):
                vertextemp[i, :] = np.array(nodeSet[sort[i]].coordinates)
            vertextemp = vertextemp[vertextemp[:, 1].argsort()]

            # split element nodes by y coordinate then sort by z coordinate then split again
            vertexint = vertextemp[0:4]
            vertexint = vertexint[vertexint[:, 2].argsort()]
            vertexlow = np.vstack((vertexlow, vertexint[0:2]))
            vertexhigh = np.vstack((vertexhigh, vertexint[2:4]))
            vertexint = vertextemp[4:8]
            vertexint = vertexint[vertexint[:, 2].argsort()]
            vertexlow = np.vstack((vertexlow, vertexint[0:2]))
            vertexhigh = np.vstack((vertexhigh, vertexint[2:4]))

            vertexlow = np.delete(vertexlow, 0, 0)
            vertexhigh = np.delete(vertexhigh, 0, 0)

            # sort in y again
            vertexlow = vertexlow[vertexlow[:, 1].argsort()]
            vertexhigh = vertexhigh[vertexhigh[:, 1].argsort()]

            # create rotational sorting in x
            vertexlow1 = vertexlow[0:2, :]
            vertexlow2 = vertexlow[2:4, :]
            vertexlow2 = vertexlow2[vertexlow2[:, 0].argsort()]
            vertexlow1 = vertexlow1[vertexlow1[:, 0].argsort()[::-1]]
            vertexlow = np.vstack((vertexlow1, vertexlow2))
            vertexhigh1 = vertexhigh[0:2, :]
            vertexhigh2 = vertexhigh[2:4, :]
            vertexhigh2 = vertexhigh2[vertexhigh2[:, 0].argsort()]
            vertexhigh1 = vertexhigh1[vertexhigh1[:, 0].argsort()[::-1]]
            vertexhigh = np.vstack((vertexhigh1, vertexhigh2))

            vertexsubply = np.vstack((vertexlow, vertexhigh))
        # sometimes Abaqus outputs nodes with 0x, y or z coordinate as something e-15
        # i.e. very close but not quite zero, this causes the code to break so set these
        # values to zero
        for i in range(0, len(vertexsubply)):
            for j in range(0, 3):
                if vertexsubply[i, j] < tol and vertexsubply[i, j] > -tol:
                    vertexsubply[i, j] = 0
        # find location of each node with respect to the cut plane
        d = np.zeros((8, 1))
        for i in range(0, 8):
            d[i] = -np.dot(perp, vertexsubply[i, :]) + np.dot(perp, r)
        # set variable value to latter check the cut plane intersects the model
        for i in range(0, len(d)):
            if d[i] > 0:
                sympos = 1
            elif d[i] < 0:
                symneg = 1
        # iterate through each node to see if any lie on the cut plane
        l = 0
        for j in range(0, 8):
            if d[j] < 0 or d[j] > 0:
                l = l + 1
            else:
                l = l
        if l != 0:
            intersection = 0
        recordj1 = []
        for j1 in range(0, 8):
            # if node lies on cut plane assign to point
            if d[j1] == 0:
                intersection += 1
                point[intersection - 1, :] = vertexsubply[j1, :]
                recordj1.append(j1)
            # if two nodes lie either side of the cut plane find the intersection point of the edge
            # between the two nodes and assign to point
            else:
                if j1 == 0:
                    for j2 in (1, 3, 4):
                        if d[j1] * d[j2] < 0:
                            intersection += 1
                            point[intersection - 1, :] = (d[j1] * vertexsubply[j2, :] - d[j2] * vertexsubply[j1, :]) / (
                                        d[j1] - d[j2])
                if j1 == 1:
                    for j2 in (2, 5):
                        if d[j1] * d[j2] < 0:
                            intersection += 1
                            point[intersection - 1, :] = (d[j1] * vertexsubply[j2, :] - d[j2] * vertexsubply[j1, :]) / (
                                        d[j1] - d[j2])
                if j1 == 2:
                    for j2 in (3, 6):
                        if d[j1] * d[j2] < 0:
                            intersection += 1
                            point[intersection - 1, :] = (d[j1] * vertexsubply[j2, :] - d[j2] * vertexsubply[j1, :]) / (
                                        d[j1] - d[j2])
                if j1 == 3:
                    j2 = 7
                    if d[j1] * d[j2] < 0:
                        intersection += 1
                        point[intersection - 1, :] = (d[j1] * vertexsubply[j2, :] - d[j2] * vertexsubply[j1, :]) / (
                                    d[j1] - d[j2])
                if j1 == 4:
                    for j2 in (5, 7):
                        if d[j1] * d[j2] < 0:
                            intersection += 1
                            point[intersection - 1, :] = (d[j1] * vertexsubply[j2, :] - d[j2] * vertexsubply[j1, :]) / (
                                        d[j1] - d[j2])
                if j1 == 5:
                    j2 = 6
                    if d[j1] * d[j2] < 0:
                        intersection += 1
                        point[intersection - 1, :] = (d[j1] * vertexsubply[j2, :] - d[j2] * vertexsubply[j1, :]) / (
                                    d[j1] - d[j2])
                if j1 == 6:
                    j2 = 7
                    if d[j1] * d[j2] < 0:
                        intersection += 1
                        point[intersection - 1, :] = (d[j1] * vertexsubply[j2, :] - d[j2] * vertexsubply[j1, :]) / (
                                    d[j1] - d[j2])
        area = 0
        angles = {}
        point = np.delete(point, np.s_[intersection:12], 0)
        pointsorted = np.zeros((intersection, 3))
        # calculate area of triangle if 3 intersection points
        if intersection == 3:
            pointsorted = point.copy()
            area = 0.5 * np.linalg.norm(np.cross((point[1, :] - point[0, :]), (point[2, :] - point[0, :])))
        # sort points if there are more than 3 intersection points
        elif intersection >= 4:
            # define reference vector
            ref = point[1, :] - point[0, :]
            # set the reference angle of the reference vector to zero
            angles[1] = 0
            for i in range(2, intersection):
                # find all other vectors with respect to the first point
                # find their angles with respect to the reference angle
                # then arrange from smallest to largest angle
                ref2 = point[i, :] - point[0, :]
                cros = np.cross(ref, ref2)
                # refNorm = np.linalg.norm(cros)
                product = np.linalg.norm(ref) * np.linalg.norm(ref2)
                dot = np.dot(ref, ref2)
                cos = dot / product

                if abs(point[0, dir] - point[1, dir]) < tol and point[i, dir] > point[0, dir]:
                    angles[i] = math.acos(cos)
                elif abs(point[0, dir] - point[1, dir]) < tol and point[i, dir] < point[0, dir]:
                    angles[i] = -math.acos(cos)
                elif cros[2] < 0:
                    angles[i] = -math.acos(cos)
                else:
                    angles[i] = math.acos(cos)
            # create sorted intersection points array
            sorted_angles = sorted(angles.items(), key=lambda x: x[1])
            sorted_index = np.array(sorted_angles)[:, 0]
            pointsorted[0, :] = point[0, :].copy()
            for i in range(0, len(sorted_index)):
                pointsorted[i + 1, :] = point[sorted_index[i], :]
        # find area of polygon if 5 or more intersection points
        if intersection >= 5:
            area = poly5(pointsorted, tol, dir)
        # find area of quadrilateral if 4 intersection points
        if intersection == 4:
            area, point1 = poly4(point)
        # if more than 3 intersection points find how plies intersect the cut plane area
        if intersection >= 3:
            t1 = 0
            for z in range(0, len(Vert)):
                c2 = 0
                pointarea = 0
                # if all intersection points lie in ply
                if np.all(pointsorted[:, dir] >= Vert[z - 1] - tol) and np.all(
                        point[:, dir] <= Vert[z] + tol) and z != 0:
                    area2[label].append(area)
                    continue
                # if all intersection points lie in ply
                elif np.all(pointsorted[:, dir] >= minVert - tol) and np.all(point[:, dir] <= Vert[z] + tol) and z == 0:
                    area2[label].append(area)
                    continue
                # if the intersection points lie outside the ply set the area to zero for that ply
                elif np.all(pointsorted[:, dir] > Vert[z] - tol):
                    area2[label].append(0)
                # if the intersection points lie outside the ply set the area to zero for that ply
                elif np.all(pointsorted[:, dir] < Vert[z - 1] + tol) and z != 0:
                    area2[label].append(0)
                elif np.all(pointsorted[:, dir] <= Vert[z] + tol) and z != 0:
                    pointtemp = np.zeros((1, 3))
                    for i in range(0, len(pointsorted)):
                        if pointsorted[i, dir] > Vert[z - 1] + tol:
                            pointtemp = np.vstack((pointtemp, pointsorted[i, :]))
                    pointtemp = np.delete(pointtemp, 0, 0)
                    pointarea = np.vstack((interceptprev, pointtemp))
                    c2 = 1
                # if ply intersects previously calculated cut plane area
                else:
                    c2 = 1
                    # initialise variables
                    intercept = np.zeros((2, 3))
                    intercept_index = []
                    d1 = np.zeros((intersection, 1))

                    for i in range(0, intersection):
                        d1[i] = -np.dot(perp1[z, :], pointsorted[i, :]) + np.dot(perp1[z, :], x[z, :])
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
                                intercept[intersection1 - 1, :] = pointsorted[j1, :]
                            else:
                                for j2 in range(j1 + 1, intersection):
                                    if d1[j1] * d1[j2] < 0:
                                        intercept_index.append(j1)
                                        intercept_index.append(j2)
                                        intersection1 += 1
                                        intercept[intersection1 - 1, :] = (d1[j1] * pointsorted[j2, :] - d1[
                                            j2] * pointsorted[
                                                  j1,
                                                  :]) / (
                                                                                  d1[j1] - d1[j2])
                    # if there are four intersection points enter this loop to find the points of intersection
                    # between the ply top surface and the cut plane area within the element
                    if intersection >= 4:
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
                                intercept[intersection1 - 1, :] = pointsorted[j1, :]
                            else:
                                if d1[j1] * d1[j1 + 1] < 0:
                                    intersection1 += 1
                                    intercept_index.append(j1)
                                    intercept_index.append(j1 + 1)
                                    intercept[intersection1 - 1, :] = (d1[j1] * pointsorted[j1 + 1, :] - d1[
                                        j1 + 1] * pointsorted[j1, :]) / (d1[j1] - d1[j1 + 1])
                        # loop round to compare last and 1st intersection points between the cut plane and element edges
                        if d1[intersection - 1] == 0:
                            intersection1 += 1
                            intercept[intersection1 - 1, :] = pointsorted[intersection - 1, :]
                        else:
                            if d1[intersection - 1] * d1[0] < 0:
                                intersection1 += 1
                                intercept_index.append(intersection - 1)
                                intercept_index.append(0)
                                intercept[intersection1 - 1, :] = (d1[intersection - 1] * pointsorted[0, :] - d1[
                                    0] * pointsorted[
                                         intersection - 1,
                                         :]) / (
                                                                          d1[intersection - 1] - d1[0])
                    # find points making up the cut plane area if no cut plane area has been calculated yet
                    if t1 == 0:
                        pointbelow = np.zeros((1, 3))
                        for i in range(0, intersection):
                            if pointsorted[i, dir] < Vert[z]:
                                pointbelow = np.vstack((pointbelow, pointsorted[i, :]))
                        pointbelow = np.delete(pointbelow, 0, 0)

                        pointarea = np.vstack((pointbelow, intercept))
                        t1 += 1
                    # if some cut plane area has been found find points making up cut plane area between the
                    # current ply plane and previous ply plane
                    else:
                        c1 = 0
                        pointinter = np.zeros((1, 3))
                        for i in range(0, intersection):
                            if pointsorted[i, dir] < Vert[z] and pointsorted[i, dir] > Vert[z - 1]:
                                pointinter = np.vstack((pointinter, pointsorted[i, :]))
                                c1 += 1
                        pointinter = np.delete(pointinter, 0, 0)
                        if c1 == 0:
                            pointarea = np.vstack((interceptprev, intercept))
                        elif c1 > 0:
                            pointarea = np.vstack((interceptprev, pointinter))
                            pointarea = np.vstack((pointarea, intercept))
                    interceptprev = intercept.copy()
                # if cut plane area has been cut by plies calculate are of triangle or
                # polygon with 4 or more nodes
                if c2 == 1 and len(pointarea) == 3:
                    area = 0.5 * np.linalg.norm(
                        np.cross((pointarea[1, :] - pointarea[0, :]), (pointarea[2, :] - pointarea[0, :])))
                    area2[label].append(area)
                elif c2 == 1 and len(pointarea) >= 4:
                    area = poly5(pointarea, tol, dir)
                    area2[label].append(area)
        # if no area exists set to 0
        if area2[label] == []:
            area2[label] = [0] * nPlies
        # if cut plane lies on the face of the element, halve the cut plane area in the element
        # to ensure the overall cut plane area is not twice what it should be
        if recordj1 == [0, 1, 2, 3] or recordj1 == [4, 5, 6, 7] or recordj1 == [0, 1, 4, 5] or recordj1 == [1, 2, 5, 6] or recordj1 == [2, 3, 6, 7] or recordj1 == [0, 3, 4, 7]:
            for j in range(0, nPlies):
                area2[label][j] = 0.5 * area2[label][j]
    else:
        area2[label] = [0] * nPlies
    area2[label].reverse()
    return(area2, sympos, symneg)