# function to calculate the cut plane area within a model meshed with tetrahedral and hexahedral elements
# import functions used within this function
from Plycomptetel import plycomptetel
from Plycomphexel import plycomphexel
import numpy as np
import math

def plycomphextet(nodeSet, p, perp, r, tol, area2, d) :
    # initialising variables
    area3 = {}
    areaelply = {}
    m = 0
    m1 = d.copy()
    t = r
    vertexsubply = np.zeros((8, 3))
    nLayups = len(p.compositeLayups.keys())
    sympos = 0
    symneg = 0

    # iterate through all the layups
    for n in range(0, nLayups):
        key = p.compositeLayups.keys()[n]
        # import number of plies in layup
        nPlies = len(p.compositeLayups[key].plies)
        plyRegion = p.compositeLayups[key].plies[0].region[0]
        # import element information
        plyElements = p.sets[plyRegion].elements

        # ply stack direction
        stackAxis = p.compositeLayups[key].orientation.stackDirection
        plyThickness = 0
        region = []

        # iterate through plies to check thicknesses add to 1
        for j in range(0, nPlies):
            plyThickness += p.compositeLayups[key].plies[j].thickness
            region.append(p.compositeLayups[key].plies[j].region[0])
            # check if relative thickness add to one
        if 0.9 > int(plyThickness) > 1.1:
            raise Exception('Defined relative thicknesses do not add up to one')

        # work out the location of ply surfaces
        Vert = np.zeros((nPlies + 1, 1))
        for j in range(0, nPlies + 1):
            for i in range(0, j):
                Vert[j] = Vert[j - 1] + p.compositeLayups[key].plies[i].thickness
        Vert = Vert[1:nPlies + 1]

        # initialise variables
        x = np.zeros((len(Vert), 3))
        perp1 = np.zeros((len(Vert), 3))
        Vertices = np.zeros((len(nodeSet), 1))
        Vertices1 = np.zeros((len(nodeSet), 1))
        Vertices2 = np.zeros((len(nodeSet), 1))

        # assign x, y and z coordinates of each vertex to variables
        for i in range(0, len(p.nodes)):
            Vertices[i] = nodeSet[i].coordinates[0]
            Vertices1[i] = nodeSet[i].coordinates[1]
            Vertices2[i] = nodeSet[i].coordinates[2]

        # if the stacking direction is in the x direction
        if str(stackAxis) == 'STACK_1':
            # find ply surface x coordinates
            Vert = (max(Vertices) - min(Vertices)) * Vert + min(Vertices)
            minVert = min(Vertices)
            # set relevant counters for later in the code
            dir = 0
            r = 0
            f = 2
            g = r
            h = 1
            # iterate through the ply surfaces to find the normal vector and point on the plane for each surface
            for i in range(0, len(Vert)):
                u = np.array([0, np.max(Vertices1), np.max(Vertices2)])
                v = np.array([0, np.min(Vertices1), np.min(Vertices2)])
                w = np.array([0, np.max(Vertices1) * 0.26, np.max(Vertices2) * 0.75])
                # np.put because Vert[i] is a list
                np.put(u, [0], Vert[i])
                np.put(v, [0], Vert[i])
                np.put(w, [0], Vert[i])
                x[i, :] = w
                pq = u - v
                pr = w - v
                perp1[i, :] = np.cross(pq, pr)

        # if the stacking direction is in the y direction
        elif str(stackAxis) == 'STACK_2':
            # find ply surface y coordinates
            Vert = (max(Vertices1) - min(Vertices1)) * Vert + min(Vertices1)
            minVert = min(Vertices1)
            # set relevant counters for later in the code
            dir = 1
            r = 1
            f = 0
            g = r
            h = 2
            # iterate through the ply surfaces to find the normal vector and point on the plane for each surface
            for i in range(0, len(Vert)):
                u = np.array([np.max(Vertices), Vert[i], np.max(Vertices2)])
                v = np.array([np.min(Vertices), Vert[i], np.min(Vertices2)])
                w = np.array([np.max(Vertices) * 0.26, Vert[i], np.max(Vertices2) * 0.75])
                x[i, :] = w
                pq = u - v
                pr = w - v
                perp1[i, :] = np.cross(pq, pr)

        # if the stacking direction is in the z direction
        elif str(stackAxis) == 'STACK_3':
            Vert = (max(Vertices2) - min(Vertices2)) * Vert + min(Vertices2)
            minVert = min(Vertices2)
            # set relevant counters for later in the code
            dir = 2
            r = 2
            f = 1
            g = r
            h = 0
            # iterate through the ply surfaces to find the normal vector and point on the plane for each surface
            for i in range(0, len(Vert)):
                u = np.array([np.max(Vertices), np.max(Vertices1), Vert[i]])
                v = np.array([np.min(Vertices), np.min(Vertices1), Vert[i]])
                w = np.array([np.max(Vertices) * 0.26, np.max(Vertices1) * 0.75, Vert[i]])
                x[i, :] = w
                pq = u - v
                pr = w - v
                perp1[i, :] = np.cross(pq, pr)

        # iterate through elements in layup
        for e in range(0, len(plyElements)) :
            # find element type
            elementType = str(p.elements[e].type)
            # run plycomphexel if element is hexahedral
            if elementType == 'C3D8R' or elementType == 'C3D20R' :
                area2, sympos, symneg = plycomphexel(e, nodeSet, perp, area2, plyElements, stackAxis, nPlies, m1, minVert, Vert, tol, vertexsubply, perp1, x, dir, t, sympos, symneg)
                areaelply = area2.copy()
            # run plycomptetel if element is tetrahedral
            if elementType == 'C3D4' or elementType == 'C3D10':
                areaelply, sympos, symneg = plycomptetel(e, nodeSet, p, perp, t, tol, area3, areaelply, Vert, minVert, r, f, g, h, x,
                             perp1, nPlies, m, sympos, symneg, m1)
                area2 = areaelply.copy()
        Totalarea = [0] * nPlies
        # calculate the total area of the cut plane in the model
        for i in range(0, nPlies):
            for j in range(0, len(p.elements)):
                Totalarea[i] = Totalarea[i] + areaelply[j + 1][i]
        # if no area is found then the plane doesn't bisect the model
        if Totalarea == [0] * nPlies :
            raise Exception('Cut plane does not bisect any of the model')
        # if any of the elements have no area in them delete the area variable associated with that element
        for label in range(0, len(p.elements)):
            if areaelply[label + 1] == [0] * nPlies:
                del areaelply[label + 1]
    # raise exception if cut plane doesn't cut the model
    if sympos == 0 or symneg == 0:
        raise Exception('Cut plane does not bisect any of the model')
    return (key, areaelply, Totalarea)