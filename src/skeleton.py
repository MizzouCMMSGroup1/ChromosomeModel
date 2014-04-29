# some imports we use
import numpy
import random
import math

# our data dimension
NUMBER_CONTACTS = 157

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

coordinate_data = []

def init_model(bounding_box=0.5):
    global coordinate_data
    for i in range(0,NUMBER_CONTACTS):
        coordinate_data += [(bounding_box*random.random(), bounding_box*random.random(), bounding_box*random.random())]


def print_model():
    global coordinate_data
    print coordinate_data


def max_if(i,j):
    return max(if_data[i,j], if_data[j,i])

def distance_sq(i,j):
    a = coordinate_data[i]
    b = coordinate_data[j]
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


def main():
    
    init_model()
    print model_score()
    #print_model()
    
    
    # SA
    # generate initial random model
    # for test in number_tests:
        # generate neighbor
        # score two models
        # if new model is better, keep
        # otherwise, accept new model based on temperature
    # output best model

    # hill climbing
    # MCMC: ??
    

if __name__ == '__main__':
    main()
