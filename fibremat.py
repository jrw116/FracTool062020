# function that uses material properties and previous area calculations to find the fracture energy and force of a particular plane
# import functions used within this function
from abaqus import *
from abaqusConstants import *
import numpy as np
import math
import os
import csv

def fibremat(area2, p, perp, multiPlane, multiPlaneCoords, singlePlane, selectAxis, axis, axisLength, z, G0, G90, sigma1, sigma2, tau12, ax, F, U, pos, pt2, plotU, plotF, q, r, s, unit, time, t1, directory, tol, Totalarea):
    # find the number of types of layup present
    nLayups = len(p.compositeLayups.keys())
    # only enter this loop if the cut plane actually cuts the model
    if area2 != {}:
        # extracting material orientation
        Name = 'csys-plane'
        if Name in p.features.keys():
            del p.features[Name]
        if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
            # creating coordinate system for plane
            p.DatumCsysByThreePoints(origin=s, name=Name, coordSysType=CARTESIAN, point1=q, point2=s + perp)
            # extracting coordinate axis
            planeAxis2 = p.datums[p.features[Name].id].axis2
        elif singlePlane == True:
            # creating coordinate system for plane
            p.DatumCsysByThreePoints(origin=s, name=Name, coordSysType=CARTESIAN, point1=q, point2=s + perp)
            # extracting coordinate axis
            planeAxis2 = p.datums[p.features[Name].id].axis2
        materialAngle = {}
        # for each defined ply layup
        for i in range(0, nLayups):
            key = p.compositeLayups.keys()[i]
            csys = p.compositeLayups[key].orientation.localCsys
            refAxis = str(p.compositeLayups[key].orientation.axis)
            addAngle = p.compositeLayups[key].orientation.angle
            # extract layup axis
            plyaxis1 = p.datums[csys].axis1
            nPlies = len(p.compositeLayups[key].plies)
            # for every ply extract ply orientation
            thetarecord = []
            for j in range(0, nPlies):
                plyOrientation = p.compositeLayups[key].plies[j].orientation
                plyAngle = p.compositeLayups[key].plies[j].orientationValue
                plyOrient = str(p.compositeLayups[key].plies[j].axis)
                plyOrientType = str(p.compositeLayups[key].plies[j].orientationType)
                # if no material orientation is defined for each individual ply then
                # orientation will be that of the layup
                if plyOrientation == None:
                    direction1 = plyaxis1.direction
                    totalAngle = addAngle + plyAngle
                    if plyOrientType == 'ANGLE_45':
                        totalAngle = addAngle + 45
                    elif plyOrientType == 'ANGLE_90':
                        totalAngle = addAngle + 90
                    elif plyOrientType == 'ANGLE_NEG45':
                        totalAngle = addAngle - 45
                    if totalAngle != 0:
                        # direction 1 is the direction of fibre
                        if refAxis == 'AXIS_2':
                            Ry = np.array([[math.cos(totalAngle * math.pi / 180), 0, -math.sin(totalAngle * math.pi / 180)],
                                           [0, 1, 0],
                                           [math.sin(totalAngle * math.pi / 180), 0, math.cos(45 * math.pi / 180)]])
                            direction1 = np.dot(Ry, plyaxis1.direction)
                        elif refAxis == 'AXIS_3':
                            Rz = np.array([[math.cos(totalAngle * math.pi / 180), -math.sin(totalAngle * math.pi / 180), 0],
                                           [-math.sin(totalAngle * math.pi / 180), math.cos(totalAngle * math.pi / 180), 0],
                                           [0, 0, 1]])
                            direction1 = np.dot(Rz, plyaxis1.direction)
                else:
                    # if there is a defined coordinate system for each ply then change the axis defined
                    plyaxis1 = p.datums[plyOrientation].axis1
                    direction1 = plyaxis1.direction
                    totalAngle = plyAngle
                    if plyOrientType == 'ANGLE_45':
                        totalAngle = 45
                    elif plyOrientType == 'ANGLE_90':
                        totalAngle = 90
                    elif plyOrientType == 'ANGLE_NEG45':
                        totalAngle = -45
                    if plyOrient == 'AXIS_2':
                        Ry = np.array(
                            [[math.cos(totalAngle * math.pi / 180), 0, -math.sin(totalAngle * math.pi / 180)], [0, 1, 0],
                             [math.sin(totalAngle * math.pi / 180), 0, math.cos(45 * math.pi / 180)]])
                        direction1 = np.dot(Ry, plyaxis1.direction)
                    elif plyOrient == 'AXIS_3':
                        Rz = np.array([[math.cos(totalAngle * math.pi / 180), -math.sin(totalAngle * math.pi / 180), 0],
                                       [-math.sin(totalAngle * math.pi / 180), math.cos(totalAngle * math.pi / 180), 0],
                                       [0, 0, 1]])
                        direction1 = np.dot(Rz, plyaxis1.direction)
                # finding angle between direction of plane axis and fibre axis
                if multiPlane == True and selectAxis == 'Create custom axis':
                    doty = np.dot(direction1, perp)
                    producty = np.linalg.norm(direction1) * np.linalg.norm(perp)
                elif multiPlaneCoords == True :
                    doty = np.dot(direction1, perp)
                    producty = np.linalg.norm(direction1) * np.linalg.norm(perp)
                elif multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
                    doty = np.dot(direction1, planeAxis2.direction)
                    producty = np.linalg.norm(direction1) * np.linalg.norm(planeAxis2.direction)
                elif singlePlane == True:
                    doty = np.dot(direction1, planeAxis2.direction)
                    producty = np.linalg.norm(direction1) * np.linalg.norm(planeAxis2.direction)
                cosy = doty / producty
                # rounding errors cause cosy to sometimes be >1 so if cosy value is close to 1 set to 1
                if cosy < 1+tol and cosy > 1-tol:
                    theta = 0
                else :
                    theta = (180/ math.pi) * math.acos(cosy)
                if theta > 90:
                    theta = 180 - theta
                # assigning angles to angle list
                thetarecord.append(theta)
            # assigning angles to elements
            for k in range(0, len(area2.keys())):
                materialAngle[area2.keys()[k]] = thetarecord

        # calculating toughness of composite in relation to fibre orientation
        # Then storing the toughness of composite in a database if existing
        # orientation already exist
        materialG = {}
        totalU = 0
        materialSigma1 = {}
        totalF = 0
        if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
            axis.append(z)
        elif multiPlane == True and selectAxis == 'Create custom axis':
            axis.append(z * axisLength)
        elif multiPlaneCoords == True :
            axis.append(z)
        # loop through every element to find the energy and force required for each element to fail
        for element in area2.keys():
            if area2[element] != [0]*nPlies :
                sigmaMatrix = []
                theta = materialAngle[element]
                # working out toughness and strength using models
                materialG[element] = G0 * np.cos(np.multiply(theta, np.pi / 180)) + G90 * np.sin(
                    np.multiply(theta, np.pi / 180))
                sigmaFibre = np.divide(sigma1, np.power(np.cos(np.multiply(theta, np.pi / 180)), 2))
                fMatrix = np.sqrt(np.power(np.divide(np.power(np.sin(np.multiply(theta, np.pi / 180)), 2), sigma2), 2)
                                  + np.power(np.divide(
                    np.multiply(np.sin(np.multiply(theta, np.pi / 180)), np.cos(np.multiply(theta, np.pi / 180))),
                    tau12), 2))
                for strength in fMatrix:
                    if strength != 0:
                        sigmaMatrix.append(1 / strength)
                    else:
                        sigmaMatrix.append(1000000000000)
                materialSigma1[element] = np.minimum(sigmaFibre, sigmaMatrix)
                # calculate energy dissipation
                totalU += np.dot(materialG[element], area2[element])
                totalF += np.dot(materialSigma1[element], area2[element])
        # if multiple cut planes are defined
        if multiPlane == True  or multiPlaneCoords == True:
            # total up the energy and force calculations
            if totalU not in U.keys():
                U[totalU] = [p.features.keys()[-2]]
                F[totalU] = [totalF]
                pos[totalU] = [z]
            else:
                U[totalU].append(p.features.keys()[-2])
                F[totalU].append(totalF)
                pos[totalU].append(z)
            # print the appropriate variables depending on the type of multiplane chosen
            if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
                print('At {0} = {1}'.format(ax, z))
                print('Total energy dissipated is: {0} J'.format(format(float(totalU), 'E')))
                print('Force required: {0} N'.format(format(float(totalF), 'E')))
                if unit == 'mm':
                    print('The cut plane area is: {0} mm^2'.format(sum(Totalarea)))
                else :
                    print('The cut plane area is: {0} m^2'.format(sum(Totalarea)))
            elif multiPlane == True and selectAxis == 'Create custom axis':
                print('At Point: {0}'.format(pt2))
                print('Total energy dissipated is: {0} J'.format(format(float(totalU), 'E')))
                print('Force required: {0} N'.format(format(float(totalF), 'E')))
                if unit == 'mm':
                    print('The cut plane area is: {0} mm^2'.format(sum(Totalarea)))
                else :
                    print('The cut plane area is: {0} m^2'.format(sum(Totalarea)))
            elif multiPlaneCoords == True :
                print('At plane number: {0}'.format(z))
                print('Total energy dissipated is: {0} J'.format(format(float(totalU), 'E')))
                print('Force required: {0} N'.format(format(float(totalF), 'E')))
                if unit == 'mm':
                    print('The cut plane area is: {0} mm^2'.format(sum(Totalarea)))
                else :
                    print('The cut plane area is: {0} m^2'.format(sum(Totalarea)))
                # output energy and force value to abaqus
            plotU.append(totalU)
            plotF.append(totalF)
        # print appropriate variables if single plane chosen and also write variables to csv file as plot function
        # doesn't run for single plane
        elif singlePlane == True:
            print('The following estimates were calculated')
            if unit == 'mm':
                print('Total energy dissipated from failure: {0} J.'.format(totalU))
                print('The cut plane area is: {0} mm^2'.format(sum(Totalarea)))
            else:
                print('Total energy dissipated from failure: {0} J.'.format(format(float(totalU),'E')))
                print('The cut plane area is: {0} m^2'.format(sum(Totalarea)))
            print('Force required: {0} N'.format(format(float(totalF),'E')))
            print('The selected fracture plane has coordinates {0},{1},{2}'.format(s,q,r))
            t2 = time.time()
            print('Run time: {0}'.format(t2-t1))
            # print failure energy, failure force, cut plane area and coordinates of cut plane to .csv file if
            # single cut plane is defined
            if unit == 'mm':
                with open(directory, 'wb') as file:
                    writer = csv.writer(file)
                    writer.writerow(["Failure energy (J)", "Failure force (N)","Cut plane area (mm^2)","Fracture plane x coordinates","Fracture plane y coordinates","Fracture plane z coordinates"])
                    writer.writerow([totalU, totalF, sum(Totalarea), s[0], s[1], s[2]])
                    writer.writerow(["", "", "", q[0], q[1], q[2]])
                    writer.writerow(["", "", "", r[0], r[1], r[2]])
            else :
                with open(directory, 'wb') as file:
                    writer = csv.writer(file)
                    writer.writerow(["Failure energy (J)", "Failure force (N)","Cut plane area (m^2)","Fracture plane x coordinates","Fracture plane y coordinates","Fracture plane z coordinates"])
                    writer.writerow([totalU, totalF, sum(Totalarea), s[0], s[1], s[2]])
                    writer.writerow(["", "", "", q[0], q[1], q[2]])
                    writer.writerow(["", "", "", r[0], r[1], r[2]])
    return (U, F, plotU, plotF, pos, axis)

