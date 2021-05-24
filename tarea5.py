'''This script calculates the equivalent
inductance matrix corresponding to
a Linnet transmission line

@Author: Alejandro Zuniga Perez'''
import numpy as np  
import math 
from scipy.constants import mu_0, epsilon_0
#Define certain parameters (magnitudes are in meters)
r_Linnet = 9.15e-03 # Data sheet
D_12 = D_13 = 12 
D_23 = 10 
d = 12e-02
GMD = np.cbrt(D_12*D_13*D_23)
mu0 = mu_0
def var_used(radius, dist):
    '''this matrix receives two parameters and defines certain variables'''
    req_p = 1.09*(math.sqrt(math.sqrt(0.778*radius*math.pow(dist,3))))
    req_p = math.log(1/req_p)
    rep = math.log(1/D_12)
    re_1 = math.log(1/D_23)
    return req_p, rep, re_1
def induc():
    '''This is the main function for calculating the inductance matrix
    for a transmission line used in file called Tarea 5'''
    req_p, rep, re_1 = var_used(r_Linnet, d)
    Mat1 = np.array([[0, rep, rep],
            [0, 0, re_1],
            [0, 0, 0]],dtype = float)
    a = np.array([req_p, req_p, req_p])
    Mat2 = np.transpose(Mat1)
    for i in range(3):
        for k in range(3):
            if i == k:
                Mat2[i,k] = 0
    Mat = Mat2 + Mat1
    for i in range(3):
        for k in range(3):
            if i == k:
                Mat[i,k] = a[i]
    Mat = (mu_0/2*np.pi)*Mat
    print('Matriz de inductancias\n',Mat)
def trasposed_inductance():
    '''This function calculates the equivalent inductance if the matrix 
    is transposed'''
    req_p = 1.09*(math.sqrt(math.sqrt(0.778*r_Linnet*math.pow(d,3))))
    v = GMD/req_p
    in_prima = (mu0/2*np.pi)*math.log(v)
    return in_prima
def capacitance(D_12, D_23):

    '''This function defines the relative heights corresponding to 
        the calculations of a transmission line projection'''
    req_p = 1.09*(math.sqrt(math.sqrt(0.778*r_Linnet*math.pow(d,3))))
    h_2 = h_3 = 2*21
    h_1 = (h_2 + 2*math.sqrt(math.pow(12,2) - math.pow(5,2)))
    h_12 = h_13 = math.sqrt((math.pow((h_2 + math.sqrt(math.pow(D_12, 2) - math.pow((D_23/2),2))),2)+math.pow((D_23/2),2)))
    h_23 = math.sqrt(math.pow(h_2,2) + math.pow(D_23, 2))
    fact_1 = math.log((np.cbrt(12*10*12))/(req_p))
    fact_2 = math.log((np.cbrt(h_2*h_3*h_1))/(np.cbrt(h_12*h_13*h_23)))
    cap = 5.55e-011*(1/(fact_1+fact_2))
    return cap 

induc()
print('Inductancia', trasposed_inductance())
print(capacitance(D_12, D_23))
