# function to sort intersection points depending on position relative other points
# import functions used within this function
import numpy as np
import math

def line(point5,point3,point,tol,i,f,g,h) :
    # find where points are positioned in relation to each other
    if abs(point3[0, f] - point3[1, f]) < tol:
        if point3[0, h] > point3[1, h]:
            maxlower = point3[0, :].copy()
            minlower = point3[1, :].copy()
        elif point3[0, h] < point3[1, h]:
            maxlower = point3[1, :].copy()
            minlower = point3[0, :].copy()

        if point3[2, h] > point3[3, h]:
            maxupper = point3[2, :].copy()
            minupper = point3[3, :].copy()
        elif point3[2, h] < point3[3, h]:
            maxupper = point3[3, :].copy()
            minupper = point3[2, :].copy()
        # find equation of max line
        if maxupper[h] == maxlower[h]:
            zm1 = maxlower[h].copy()
        else:
            # max line
            m1 = (maxupper[g] - maxlower[g]) / (maxupper[h] - maxlower[h])
            cm1 = maxupper[g] - m1 * maxupper[h]
            zm1 = (point[i, g] - cm1) / m1
        # find equation of min line
        if minupper[h] == minlower[h]:
            zm2 = minlower[h].copy()
        else:
            # min line
            m2 = (minupper[g] - minlower[g]) / (minupper[h] - minlower[h])
            cm2 = minupper[g] - m2 * minupper[h]
            zm2 = (point[i, g] - cm2) / m2
        # create array of ordered points
        if point[i, h] > zm1 and point[i, h] > zm2:
            point5[0, :] = maxlower.copy()
            point5[1, :] = maxupper.copy()
        elif point[i, h] < zm1 and point[i, h] < zm2:
            point5[0, :] = minlower.copy()
            point5[1, :] = minupper.copy()

    else:
        if point3[0, f] >= point3[1, f]:
            maxlower = point3[0, :].copy()
            minlower = point3[1, :].copy()
        elif point3[0, f] < point3[1, f]:
            maxlower = point3[1, :].copy()
            minlower = point3[0, :].copy()

        if point3[2, f] >= point3[3, f]:
            maxupper = point3[2, :].copy()
            minupper = point3[3, :].copy()
        elif point3[2, f] < point3[3, f]:
            maxupper = point3[3, :].copy()
            minupper = point3[2, :].copy()
        # find equation of max line
        if maxupper[f] == maxlower[f]:
            xm1 = maxlower[f].copy()
        else:
            # max line
            m1 = (maxupper[g] - maxlower[g]) / (maxupper[f] - maxlower[f])
            cm1 = maxupper[g] - m1 * maxupper[f]
            xm1 = (point[i, g] - cm1) / m1
        # find equation of min line
        if minupper[f] == minlower[f]:
            xm2 = minlower[f].copy()
        else:
            # min line
            m2 = (minupper[g] - minlower[g]) / (minupper[f] - minlower[f])
            cm2 = minupper[g] - m2 * minupper[f]
            xm2 = (point[i, g] - cm2) / m2
        # create array of ordered points
        if point[i, f] > xm2 and point[i, f] > xm2:
            point5[0, :] = maxlower.copy()
            point5[1, :] = maxupper.copy()
        elif point[i, f] < xm2 and point[i, f] < xm2:
            point5[0, :] = minlower.copy()
            point5[1, :] = minupper.copy()
    return(point5)