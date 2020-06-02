"""
Created on Tuesday June 2nd 2020

Authors: Jack Wheaton, Johnson Lee
"""

from abaqus import *
from abaqusConstants import *
from axisplanedel import axisplanedel
from plycomp import plycomp
from fibremat import fibremat
from plot import plot
from plycomptet import plycomptet
from plycomphextet import plycomphextet
from plycomp import plycomp
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import numpy as np
import math
import time
from collections import Counter
import sys
import csv
import os

def fEDC(part, model, singlePlane, multiPlane, multiPlaneCoords, **kwargs):
    # import variables from user interface as defined in plugin3DB.py
    unit = kwargs.get('unit')
    mat = kwargs.get('mat')
    cusMat = kwargs.get('cusMat')
    N = kwargs.get('N')
    g0 = kwargs.get('g0')
    g90 = kwargs.get('g90')
    cusSigma = kwargs.get('cusSigma')
    cusSigma2 = kwargs.get('cusSigma2')
    cusTau = kwargs.get('cusTau')
    tol = kwargs.get('tolerance')
    selectAxis = kwargs.get('selectAxis')
    axisGlobal = kwargs.get('axisGlobal')
    directory1 = kwargs.get('directory')
    fileName = kwargs.get('fileName')

    #Signal start of new simulation
    print('-----------------------------New Simulation-----------------------------')

    # import model and part information
    p = mdb.models[model].parts[part]

    # getting nodes, elements and element types from part
    nodeSet = p.nodes
    elementSet = p.elements

    # run the apdelete function in order to delete any old planes or axes in the viewport
    axisplanedel(model, part)

    # ensure the correct directory is set for the Results.csv file to be saved to
    if directory1 == 'FracTool Plug-in Directory':
        dirname = os.path.dirname(os.path.abspath(__file__))
        directory = os.path.join(dirname, 'Results.csv')
    else :
        directory = os.path.join(directory1, 'Results.csv')

    # exception handling
    # check if part exists
    if part == '':
        raise Exception('Part not defined please make a part before using software')
    # check if input for number of fracture planes exist
    if N == None:
        raise Exception('Enter number of fracture planes')
    if g0 == None or g90 == None:
        raise Exception('Enter properties of material')
    if len(p.compositeLayups.keys()) == 0:
        raise Exception('Please define layup before using tool')

    # check if the tab check boxes are ticked and raise an exception if not
    if multiPlane == False and singlePlane == False and multiPlaneCoords == False:
        raise Exception('No checkboxes were ticked please tick at least one checkbox')
    #check that only 1 selection tab is selected
    if multiPlane == True and singlePlane == True:
        raise Exception('More than one tab selected, only one tab should be selected')
    if multiPlane == True and multiPlaneCoords == True :
        raise Exception('More than one tab selected, only one tab should be selected')
    if singlePlane == True and multiPlaneCoords == True :
        raise Exception('More than one tab selected, only one tab should be selected')

    # timing script
    t1 = time.time()
    # initialising variables
    area2 = {}
    U = {}
    F ={}
    pos = {}
    axis = []
    plotU = []
    plotF = []
    plotArea = []
    x = []
    y = []
    z = []

    # checking that the element types used in the model meshing are supported by the plug-in
    elementTypesupport = 0
    elementTypesupport1 = 0
    elementTypesupport2 = 0
    for i in range(0,len(elementSet)) :
        elementType = p.elements[i].type
        # set counters to ensure the correct functions are run to calculate the cut plane area
        if elementType == C3D8R or elementType == C3D20R :
            elementTypesupport = 1
        elif elementType == C3D4 or elementType == C3D10 :
            elementTypesupport1 = 1
        else :
            elementTypesupport2 = 1
        if elementTypesupport2 == 1:
            raise Exception('The mesh contains element types that are not supported')
    # if no nodes are present give error
    if len(elementSet) == 0 or len(elementSet) == None:
        raise Exception('No mesh defined please define mesh before using software')

    # setting material properties from database
    if cusMat == 'Material from database':
        if mat == 'T800s/M21':
            G90 = 0.255
            G0 = 209
            sigma1 = 2451
            sigma2 = 147
            tau12 = 145
        elif mat == 'T300/913':
            G90 = 0.211
            G0 = 91.6
            sigma1 = 582
            sigma2 = 65.5
            tau12 = 67.9
        elif mat == 'T300/920':
            G90 = 0.456
            G0 = 132
            sigma1 = 582
            sigma2 = 34.9
            tau12 = 67.9
        elif mat == 'TR50s/K51' :
            G90 = 0.255
            G0 = 32.2
            sigma1 = 2500
            sigma2 = 60
            tau12 = 88.5
        elif mat == 'IM7/8551-7' :
            G90 = 0.450
            G0 = 112.7 #TBC
            sigma1 = 2760
            sigma2 = 75.8
            tau12 = 100
    # custom material property setting
    if cusMat == 'Custom Material':
        G90 = g90
        G0 = g0
        sigma1 = cusSigma
        sigma2 = cusSigma2
        tau12 = cusTau
    # converting units if geometry of model is in millimetres
    if unit == 'mm':
        G90 = G90/1000
        G0 = G0/1000
    # converting units if geometry of model is in metres
    elif unit == 'm':
        G90 = G90*1000
        G0 = G0*1000
        sigma1 = sigma1*1000000
        sigma2 = sigma2*1000000
        tau12 = tau12*1000000

    # when multiple cut planes along the x,y or z axis are being used
    if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
        # getting coordinates of the model
        for n in range(0,len(nodeSet)):
            x.append(nodeSet[n].coordinates[0])
            y.append(nodeSet[n].coordinates[1])
            z.append(nodeSet[n].coordinates[2])
        # finding maximum and minimum coordinates of the model
        minX = min(x)
        minY = min(y)
        minZ = min(z)
        maxX = max(x)
        maxY = max(y)
        maxZ = max(z)
        # setting the start and end values for the planes
        if axisGlobal == 'X':
            start = minX+1
            end = maxX+1
        elif axisGlobal == 'Y':
            start = minY+1
            end = maxY+1
        elif axisGlobal == 'Z': 
            start = minZ+1
            end = maxZ+1
        #iterating through planes
        for z in np.arange(start,end,abs(start-end)/N):
            # Signal start of new plane
            print('-----------------------------New Plane-----------------------------')
            # Set the 3 coordinates that define the cut plane
            if axisGlobal == 'X':
                s = np.array([z,0,0])
                q = np.array([z,0,1])
                r = np.array([z,1,0])
                ax = 'x'
            elif axisGlobal == 'Y':
                s = np.array([0,z,0])
                q = np.array([0,z,1])
                r = np.array([1,z,0])
                ax = 'y'
            elif axisGlobal == 'Z':
                s = np.array([0,0,z])
                q = np.array([0,1,z])
                r = np.array([1,0,z])
                ax = 'z'
            # create plane for visualisation in abaqus cae
            p.DatumPlaneByThreePoints(point1=s, 
                point2=q, 
                point3=r) 
            # calculating vectors
            pq = q - s
            pr = r - s
            # equation of plane
            perp = np.cross(pq,pr)
            d = -np.dot(perp,-s)
            # choosing the correct function to run to find the cut plane area within each element split by the amount of element area in each ply
            if elementTypesupport == 1 and elementTypesupport1== 0 :
                print('Hexahedral mesh chosen')
                area2, key, Totalarea = plycomp(p, area2, nodeSet, perp, d, tol, r)
            elif elementTypesupport == 0 and elementTypesupport1== 1 :
                print('Tetrahedral mesh chosen')
                key, areaelply, Totalarea = plycomptet(nodeSet, p, perp, r, tol, d)
                area2 = areaelply.copy()
            elif elementTypesupport == 1 and elementTypesupport1 == 1:
                print('Mixed Hex and Tet mesh')
                key, areaelply, Totalarea = plycomphextet(nodeSet, p, perp, r, tol, area2, d)
                area2 = areaelply.copy()
            pt2 = []
            axisLength = []
            plotArea.append(sum(Totalarea))
            # run the function fibremat to calculate the force and energy outputs
            U, F, plotU, plotF, pos, axis = fibremat(area2, p, perp, multiPlane, multiPlaneCoords, singlePlane, selectAxis, axis, axisLength, z, G0, G90, sigma1, sigma2, tau12, ax, F, U, pos, pt2, plotU, plotF, q, r, s, unit, time, t1, directory, tol, Totalarea)
        # run function to plot the appropriate graphs and save the data required
        plot(U, plotU, pos, F, plotF, plotArea, ax, t1, pt2, axis, multiPlane, multiPlaneCoords, selectAxis, directory, unit)

    # if coordinates of the cut planes are imported from a .csv file
    if multiPlaneCoords == True :
        # csv file import
        csv = np.genfromtxt(fileName, delimiter=',')
        s1 = csv[:, 0:3]
        q1 = csv[:, 3:6]
        r1 = csv[:, 6:9]
        #iterating through planes
        for z in range(0,len(s1)) :
            # Signal start of new plane
            print('-----------------------------New Plane-----------------------------')
            s = s1[z, :]
            q = q1[z, :]
            r = r1[z, :]
            p.DatumPlaneByThreePoints(point1=s, point2=q, point3=r)
            # calculating vectors
            pq = q - s
            pr = r - s
            # equation of plane
            perp = np.cross(pq, pr)
            d = -np.dot(perp, -s)
            # choosing the correct function to run to find the cut plane area within each element split by the amount of element area in each ply
            if elementTypesupport == 1 and elementTypesupport1== 0 :
                print('Hexahedral mesh chosen')
                area2, key, Totalarea = plycomp(p, area2, nodeSet, perp, d, tol, r)
            elif elementTypesupport == 0 and elementTypesupport1== 1 :
                print('Tetrahedral mesh chosen')
                key, areaelply, Totalarea = plycomptet(nodeSet, p, perp, r, tol, d)
                area2 = areaelply.copy()
            elif elementTypesupport == 1 and elementTypesupport1 == 1:
                print('Mixed Hex and Tet mesh')
                key, areaelply, Totalarea = plycomphextet(nodeSet, p, perp, r, tol, area2, d)
                area2 = areaelply.copy()
            pt2 = []
            axisLength = []
            ax = []
            plotArea.append(sum(Totalarea))
            # run the function fibremat to calculate the force and energy outputs
            U, F, plotU, plotF, pos, axis = fibremat(area2, p, perp, multiPlane, multiPlaneCoords, singlePlane, selectAxis, axis, axisLength, z, G0, G90, sigma1, sigma2, tau12, ax, F, U, pos, pt2, plotU, plotF, q, r, s, unit, time, t1, directory, tol, Totalarea)
        # run function to plot the appropriate graphs and save the data required
        plot(U, plotU, pos, F, plotF, plotArea, ax, t1, pt2, axis, multiPlane, multiPlaneCoords, selectAxis, directory, unit)

    # if own axis defined, the code is the same as above
    # Only difference is the plane definition
    if multiPlane == True and selectAxis == 'Create custom axis':
        axisPoints = kwargs.get('axisPoints')
        # timing script
        t1 = time.time()
        # initialising variables
        U = {}
        pos = {}
        axis = []
        plotU = []
        # selecting points n viewport
        if axisPoints == 'Select 2 points in viewport to create axis':
            axisPointStart = kwargs.get('axisPointStart') 
            axisPointEnd = kwargs.get('axisPointEnd')
            #exception handling
            if axisPointStart == None or axisPointEnd == None:
                raise Exception('Points not selected')
            s = np.array(axisPointStart.coordinates)
            q = np.array(axisPointEnd.coordinates)
        # define points manually
        elif axisPoints == 'Define points manually':
            coordAxisStart = kwargs.get('coordAxisStart')
            coordAxisEnd = kwargs.get('coordAxisEnd')
            # exception handling
            if coordAxisStart == () or coordAxisEnd == ():
                raise Exception('Points not defined')
            elif len(coordAxisStart[0]) != 3 or len(coordAxisEnd[0]) != 3:
                raise Exception('Axis points not defined properly')
            s = np.array([coordAxisStart[0][0],coordAxisStart[0][1],coordAxisStart[0][2]])
            q = np.array([coordAxisEnd[0][0],coordAxisEnd[0][1],coordAxisEnd[0][2]])
        # find the axis length and plot the normal axis to the planes in the viewport
        perp = np.subtract(q,s)
        axisLength = np.linalg.norm(perp)
        p.DatumAxisByTwoPoint(point1=s, 
            point2=q)
        d1 = p.datums
        planeName = p.features.keys()[-1]
        inc = 1/float(N)
        # iterate through the planes
        for z in np.arange(0,1, inc):
            # Signal start of new plane
            print('-----------------------------New Plane-----------------------------')
            pt2 = s + np.multiply(perp,z)
            d = np.dot(perp,pt2)
            r = pt2.copy()
            # plot plane using point on plane and normal vector to the plane
            p.DatumPlaneByPointNormal(normal=d1[p.features[planeName].id], point=s + np.multiply(perp,z))
            # choosing the correct function to run to find the cut plane area within each element split by the amount of element area in each ply
            if elementTypesupport == 1 and elementTypesupport1== 0 :
                print('Hexahedral mesh chosen')
                area2, key, Totalarea = plycomp(p, area2, nodeSet, perp, d, tol, r)
            elif elementTypesupport == 0 and elementTypesupport1 == 1:
                print('Tetrahedral mesh chosen')
                key, areaelply, Totalarea = plycomptet(nodeSet, p, perp, r, tol, d)
                area2 = areaelply.copy()
            elif elementTypesupport == 1 and elementTypesupport1 == 1:
                print('Mixed Hex and Tet mesh')
                key, areaelply, Totalarea = plycomphextet(nodeSet, p, perp, r, tol, area2, d)
                area2 = areaelply.copy()
            ax = []
            r = []
            plotArea.append(sum(Totalarea))
            # run the function fibremat to calculate the force and energy outputs
            U, F, plotU, plotF, pos, axis = fibremat(area2, p, perp, multiPlane, multiPlaneCoords, singlePlane, selectAxis, axis, axisLength, z, G0, G90, sigma1, sigma2, tau12, ax, F, U, pos, pt2, plotU, plotF, q, r, s, unit, time, t1, directory, tol, Totalarea)
        # run function to plot the appropriate graphs and save the data required
        plot(U, plotU, pos, F, plotF, plotArea, ax, t1, pt2, axis, multiPlane, multiPlaneCoords, selectAxis, directory, unit)
        
    # this if loop runs if only one cut plane is defined
    if singlePlane == True:
        selectPsingle = kwargs.get('selectPsingle')
        # select points in viewport
        if selectPsingle == 'Select points on part display':
            coord1 = kwargs.get('coord1')
            coord2 = kwargs.get('coord2')
            coord3 = kwargs.get('coord3')
            # exception handling
            if coord1 == None or coord2 == None or coord3 == None:
                raise Exception('Point not picked please pick points')
            s = np.array(coord1.coordinates)
            q = np.array(coord2.coordinates)
            r = np.array(coord3.coordinates)
        # define points manually
        elif selectPsingle == 'Define Points Manually':
            enterCoord1 = kwargs.get('enterCoord1')
            enterCoord2 = kwargs.get('enterCoord2')
            enterCoord3 = kwargs.get('enterCoord3')
            # exception handling
            if enterCoord1 == () or enterCoord2 == () or enterCoord3 == ():
                raise Exception('Points not defined please define point')
            elif len(enterCoord1[0]) != 3 or len(enterCoord2[0]) != 3 or len(enterCoord3[0]) != 3:
                raise Exception('Axis points not defined properly')
            s = np.array([enterCoord1[0][0],enterCoord1[0][1],enterCoord1[0][2]])
            q = np.array([enterCoord2[0][0],enterCoord2[0][1],enterCoord2[0][2]])
            r = np.array([enterCoord3[0][0],enterCoord3[0][1],enterCoord3[0][2]])
        # plot cut plane
        p.DatumPlaneByThreePoints(point1=s, point2=q, point3=r)
        # carry out calculations to find the normal vector to the plane
        pq = q - s
        pr = r - s
        perp = np.cross(pq,pr)
        d = -np.dot(perp,-s)
        # choosing the correct function to run to find the cut plane area within each element split by the amount of element area in each ply
        if elementTypesupport == 1 and elementTypesupport1 == 0:
            print('Hexahedral mesh chosen')
            area2, key, Totalarea = plycomp(p, area2, nodeSet, perp, d, tol, r)
        elif elementTypesupport == 0 and elementTypesupport1 == 1:
            print('Tetrahedral mesh chosen')
            key, areaelply, Totalarea = plycomptet(nodeSet, p, perp, r, tol, d)
            area2 = areaelply.copy()
        elif elementTypesupport == 1 and elementTypesupport1 == 1:
            print('Mixed Hex and Tet mesh')
            key, areaelply, Totalarea = plycomphextet(nodeSet, p, perp, r, tol, area2, d)
            area2 = areaelply.copy()
        print(Totalarea)
        #Material function
        axisLength = []
        ax = []
        pt2 = []
        # run the function fibremat to calculate the force and energy outputs and output these to a csv file
        fibremat(area2, p, perp, multiPlane, multiPlaneCoords, singlePlane, selectAxis, axis, axisLength, z, G0, G90,
                 sigma1, sigma2, tau12, ax, F, U, pos, pt2, plotU, plotF, q, r, s, unit, time, t1, directory, tol, Totalarea)
    # code end