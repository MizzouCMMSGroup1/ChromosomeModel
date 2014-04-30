# some imports we use
import numpy
import random
import math
from scipy import optimize

# our data dimension (*xyz)
NUMBER_CONTACTS = 157
NUMBER_CONTACTS_POINTS = NUMBER_CONTACTS * 3

# bring in our cleaned data
IF_FILENAME = "if_data_stripped.csv"
if_data = numpy.loadtxt(IF_FILENAME, delimiter=',')

IF_TOTAL = numpy.sum(if_data)

# weighting scores from paper, we are using megabase resolution
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
def f(x, *args):
    print x
    
    global coordinate_data
    for i in range(0,NUMBER_CONTACTS_POINTS):
        coordinate_data[i] = x[i]
    
    return model_score()

def main():
        
    random_start = init_model()
    
    #print_model()
    
    #return
    
    point_data = []
    args = []
    
    #print point_data
    
    #return
    opts = {'maxiter' : None,    # default value.
         'disp' : True,    # non-default value.
         'gtol' : 1e-5,    # default value.
         'norm' : numpy.inf,  # default value.
         'eps' : 1.4901161193847656e-08}  # default value.

    res2 = optimize.minimize(f, random_start, args=args,
                         method='CG', options=opts)
    
    print res2
    # hill climbing
    # MCMC: ??
    

if __name__ == '__main__':
    main()
