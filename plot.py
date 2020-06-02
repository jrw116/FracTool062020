# function to plot the appropriate graphs in abaqus and save data to csv files
# import functions used within this function
from abaqus import *
from abaqusConstants import *

import numpy as np
import visualization
import xyPlot
import time
import csv
import os

def plot(U, plotU, pos, F, plotF, plotArea, ax, t1, pt2, axis, multiPlane, multiPlaneCoords, selectAxis, directory, unit):
    # initialise variables
    fileNumber = 1
    critU = min(U.keys())
    planes = U[critU]
    force = F[critU]
    # if x,y or z multiplane chosen write certain values to Abaqus dialog box
    if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
        critPos = pos[critU]
        # output to message area in abaqus
        print('From the estimated values')
        print('Critical energy is: {0} .'.format(format(format(critU, 'E'))))
        print('Force required: {0} N'.format(format(force[0], 'E')))
        print('At the following planes: {0}, at {1} = {2}'.format(planes, ax, critPos))
    # if multiplane with chosen axis write certain values to Abaqus dialog box
    elif multiPlaneCoords == True :
        critPos = pos[critU]
        print('The critical fracture plane/s: {0}'.format(planes))
        print('At plane number: {0}'.format(critPos))
        print('With critical fracture energy: {0} J'.format(format(critU, 'E')))
        print('Force required: {0} N'.format(format(force[0], 'E')))
    # if multiplane with chosen axis write certain values to Abaqus dialog box
    else :
        print('The critical fracture plane/s: {0}'.format(planes))
        print('At point: {0}'.format(pt2))
        print('With critical fracture energy: {0} J'.format(format(critU, 'E')))
        print('Force required: {0} N'.format(format(force[0], 'E')))
    # print time taken to run
    t2 = time.time()
    print('Run time: {0}'.format(t2 - t1))
    # energy plot
    # plotting output energies with increasing co-ordinates
    # getting existing xyplots
    outputFiles = session.xyPlots.keys()
    fileName = 'Energy Output{0}'.format(fileNumber)
    while fileName in outputFiles:
        fileNumber += 1
        fileName = 'Energy Output{0}'.format(fileNumber)
    # creating new data plot
    xyp = session.XYPlot(fileName)
    # creating new chart
    chartName = xyp.charts.keys()[0]
    chart = xyp.charts[chartName]

    # x and y label
    yQuantity = visualization.QuantityType(type=ENERGY)
    c = np.empty((len(axis), 2))
    c[:, 0] = axis
    c[:, 1] = plotU
    # creating data plot
    xy2 = xyPlot.XYData(data=(c), sourceDescription='Entered from keyboard',
                        axis2QuantityType=yQuantity, )
    # creating curve
    c2 = session.Curve(xyData=xy2)
    # plotting curve and changing axis label
    chart.setValues(curvesToPlot=(c2,), )
    if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
        session.charts[chartName].axes1[0].axisData.setValues(useSystemTitle=False,
                                                          title='{0} coordinate'.format(ax))
    elif multiPlaneCoords == True :
        session.charts[chartName].axes1[0].axisData.setValues(useSystemTitle=False,
                                                          title='Plane Number')
    else:
        session.charts[chartName].axes1[0].axisData.setValues(useSystemTitle=False,
                                                          title='Distance along axis path')
    # changing view to graph
    session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
    # writing different variables to csv file depending on tab chosen in FracTool GUI
    # (either x,y or z or defined normal axis), .csv file also changed depending on model geometry units
    if unit == 'mm':
        if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
            with open(directory, 'wb') as file:
                writer = csv.writer(file)
                writer.writerow(["{0} coordinate".format(ax),"Cut plane area (mm^2)","Failure Energy (J)","Failure Force (N)","Critical failure energy (J)", "Critical failure force (N)"])
                writer.writerow([axis[0], plotArea[0], plotU[0], plotF[0], critU, force[0]])
                for num in range(1,len(axis)) :
                    writer.writerow([axis[num], plotArea[num], plotU[num], plotF[num], "", ""])
        elif multiPlaneCoords == True :
            with open(directory, 'wb') as file:
                writer = csv.writer(file)
                writer.writerow(["Plane Number {0}","Cut plane area (mm^2)","Failure Energy (J)","Failure Force (N)","Critical failure energy (J)", "Critical failure force (N)"])
                writer.writerow([axis[0], plotArea[0], plotU[0], plotF[0], critU, force[0]])
                for num in range(1,len(axis)) :
                    writer.writerow([axis[num], plotArea[num], plotU[num], plotF[num], "", ""])
        else :
            with open(directory, 'wb') as file:
                writer = csv.writer(file)
                writer.writerow(["Distance along axis path","Cut plane area (mm^2)","Failure Energy (J)","Failure Force (N)","Critical failure energy (J)", "Critical failure force (N)"])
                writer.writerow([axis[0], plotArea[0], plotU[0], plotF[0], critU, force[0]])
                for num in range(1,len(axis)) :
                    writer.writerow([axis[num], plotArea[num], plotU[num], plotF[num], "", ""])
    else :
        if multiPlane == True and selectAxis == 'Direction of axis perpendicular to fracture planes':
            with open(directory, 'wb') as file:
                writer = csv.writer(file)
                writer.writerow(["{0} coordinate".format(ax),"Cut plane area (m^2)","Failure Energy (J)","Failure Force (N)","Critical failure energy (J)", "Critical failure force (N)"])
                writer.writerow([axis[0], plotArea[0], plotU[0], plotF[0], critU, force[0]])
                for num in range(1,len(axis)) :
                    writer.writerow([axis[num], plotArea[num], plotU[num], plotF[num], "", ""])
        elif multiPlaneCoords == True :
            with open(directory, 'wb') as file:
                writer = csv.writer(file)
                writer.writerow(["Plane Number","Cut plane area (m^2)","Failure Energy (J)","Failure Force (N)","Critical failure energy (J)", "Critical failure force (N)"])
                writer.writerow([axis[0], plotArea[0], plotU[0], plotF[0], critU, force[0]])
                for num in range(1,len(axis)) :
                    writer.writerow([axis[num], plotArea[num], plotU[num], plotF[num], "", ""])
        else :
            with open(directory, 'wb') as file:
                writer = csv.writer(file)
                writer.writerow(["Distance along axis path","Cut plane area (m^2)","Failure Energy (J)","Failure Force (N)","Critical failure energy (J)", "Critical failure force (N)"])
                writer.writerow([axis[0], plotArea[0], plotU[0], plotF[0], critU, force[0]])
                for num in range(1,len(axis)) :
                    writer.writerow([axis[num], plotArea[num], plotU[num], plotF[num], "", ""])