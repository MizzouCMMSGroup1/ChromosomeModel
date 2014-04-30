# some imports we use
import numpy
import random
import math
from scipy import optimize
import argparse

# matplot lib
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

# our data dimension (*xyz), we are using megabase resolution
NUMBER_CONTACTS = 157
NUMBER_CONTACTS_POINTS = NUMBER_CONTACTS * 3

# bring in our cleaned data
IF_FILENAME = "if_data_stripped.csv"
if_data_raw = numpy.loadtxt(IF_FILENAME, delimiter=',')

# used to normalize our IF weights
IF_TOTAL = numpy.sum(if_data_raw)

# chromosome 7 weighting scores from paper, megabase
W1 = 1.0
W2 = 1.5
W3 = 1.5
W4 = 1.5

# micrometers
d_sq_min = 0.2
da_sq_max = 1.8
d_max = 4.5
d_sq_c = 7.0

d_sq_max = d_max * d_max # not defined in paper?

# some globals
coordinate_data = numpy.zeros(NUMBER_CONTACTS_POINTS)

def setup_testing(number_contacts=157):
    # 157 is our default/actual size of c7
    global if_data, IF_TOTAL
    global NUMBER_CONTACTS, NUMBER_CONTACTS_POINTS
    NUMBER_CONTACTS = number_contacts
    NUMBER_CONTACTS_POINTS = NUMBER_CONTACTS * 3
    if_data = if_data_raw[0:number_contacts,0:number_contacts]
    IF_TOTAL = numpy.sum(if_data)
    coordinate_data = numpy.zeros(number_contacts)


def init_model(bounding_box=0.5):
    global coordinate_data
    for i in range(0,NUMBER_CONTACTS_POINTS):
        coordinate_data[i] = bounding_box * (0.5 - random.random())
    return coordinate_data


def print_model():
    global coordinate_data
    print coordinate_data


def max_if(i,j):
    return max(if_data[i,j], if_data[j,i])


def distance_sq(i,j):
    a = [coordinate_data[i], coordinate_data[i+1], coordinate_data[i+2]]
    b = [coordinate_data[j], coordinate_data[j+1], coordinate_data[j+2]]
    return (a[0] - b[0])**2 + (a[1] - b[1])**2 + (a[2] - b[2])**2


def contact_score():
    '''
    minimize the distance (but keep above min_threshold) between non-sequential pairs that have affinity
    '''
    global IF_TOTAL
    score = 0
    for i in range(0,NUMBER_CONTACTS):
        for j in range(0,NUMBER_CONTACTS):
            if i != j or abs(i-j) != 1:
                d_sq_ij = distance_sq(i,j)
                score += W1 * math.tanh(d_sq_c - d_sq_ij) * max_if(i,j) + W2 * math.tanh(d_sq_ij - d_sq_min) / IF_TOTAL
    return score


def noncontact_score():
    '''
    maximize the distance (but keep below max_threshold) between non-sequential pairs that don't have affinity
    '''
    global IF_TOTAL
    score = 0
    for i in range(0,NUMBER_CONTACTS):
        for j in range(0,NUMBER_CONTACTS):
            if i != j or abs(i-j) != 1:
                d_sq_ij = distance_sq(i,j)
                score += W3 * math.tanh(d_sq_max - d_sq_ij) / IF_TOTAL + W4 * math.tanh(d_sq_ij - d_sq_c) / IF_TOTAL
    return score


def pair_smoothing():
    '''
    keep adjacent contacts (eg |i-j|==1) with slightly lower score than above so they are prioritized for optimization
    '''
    global IF_TOTAL
    score = 0
    for i in range(0,NUMBER_CONTACTS):
        for j in range(0,NUMBER_CONTACTS):
            if abs(i-j) == 1:
                d_sq_ij = distance_sq(i,j)
                score += W1 * max_if(i,j) * math.tanh(da_sq_max - d_sq_ij) + W2 * math.tanh(d_sq_ij - d_sq_min) / IF_TOTAL
    return score


def model_score():
    return contact_score() + noncontact_score() + pair_smoothing()


# shim between skeleton and cg code
iter_tracker = 0
old_score = 0

def f(x, *args):
    #print x
    global iter_tracker, old_score
    iter_tracker += 1
    
    global coordinate_data
    for i in range(0,NUMBER_CONTACTS_POINTS):
        coordinate_data[i] = x[i]
    
    current_score = model_score()

    #print "iter:", iter_tracker, "score:", current_score, "change:", current_score - old_score

    old_score = current_score
    
    return current_score


def main():
    global iter_tracker
    
    parser = argparse.ArgumentParser(description="Runner for PyChromosomeModeler")
    parser.add_argument('-n','--number', help='number of residues to model',type=int,default=5)
    t_group = parser.add_mutually_exclusive_group()
    t_group.add_argument('-c','--conjugate',action='store_true')
    t_group.add_argument('-a','--anneal',action='store_true')
    args = parser.parse_args()

    number_residues = args.number

    TESTING_CONGUGATE_GRADIENT = args.conjugate
    setup_testing(number_residues)

    random_start = init_model().copy()
    args = []
    
    opts = {'maxiter' : 100, 'disp' : True }

    results = 0
    
    if (TESTING_CONGUGATE_GRADIENT):
        results = optimize.minimize(f, random_start, args=args, method='CG', options=opts)    
    else:
        results = optimize.minimize(f, random_start, args=args, method='Anneal', options=opts)

    print "internal iter: ", iter_tracker

    print results
    print "saving final contact xyz coordinates"
    
    x = numpy.zeros(NUMBER_CONTACTS)
    y = numpy.zeros(NUMBER_CONTACTS)
    z = numpy.zeros(NUMBER_CONTACTS)
    
    for i in range(0,NUMBER_CONTACTS):
        x[i] = results.x[i]
        y[i] = results.x[i+1]
        z[i] = results.x[i+2]
        print results.x[i], results.x[i+1], results.x[i+2]

    mpl.rcParams['legend.fontsize'] = 10

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot(x, y, z, label='3d plot of generated contacts')
    ax.legend()

    plt.show()


if __name__ == '__main__':
    main()
