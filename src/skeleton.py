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

# some globals
coordinate_data = []
random_coordinates = []
random_offset = 0
old_coordinates = []

def init_model(bounding_box=0.5):
    global coordinate_data
    coordinate_data = []
    for i in range(0,NUMBER_CONTACTS):
        coordinate_data += [(bounding_box*random.random(), bounding_box*random.random(), bounding_box*random.random())]


def randomize_model(bounding_box=0.5):
    global coordinate_data
    global random_offset, random_coordinates, old_coordinates
    random_offset = random.randrange(0, NUMBER_CONTACTS)
    random_coordinates = (bounding_box*random.random(), bounding_box*random.random(), bounding_box*random.random())
    old_coordinates = coordinate_data[random_offset]
    coordinate_data[random_offset] = random_coordinates


def revert_model():
    print "reverting change"
    global coordinate_data
    global random_offset, old_coordinates
    coordinate_data[random_offset] = old_coordinates


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

#temperature calculator. non-linear decrease
def sigmoid_temperature(k):
  return -5000/(1 + math.exp(-k/200)) + 5000

def linear_temperature(k):
  return (-2500/1000)*k + 2500

def main():
    
    init_model()
    #print model_score()
    
    NUMBER_SIMULATIONS = 100
    TESTING_SIGMOID = True
    best_prior_score = 1e7
    
    for i in range(1,NUMBER_SIMULATIONS+1,1):
        
        T = 0
        if TESTING_SIGMOID:
            T = sigmoid_temperature(i)
        else:
            T = linear_temperature(i)
        
        randomize_model()
        
        i_score = model_score()
        print "i_score", i_score
        
        score_diff = i_score - best_prior_score
        
        if score_diff > 0:
            #print("score diff", score_diff)
            score_diff = score_diff if score_diff < 5*T else 5*T
            prob_to_accept = math.exp(-100*score_diff/T)
            print("probability to accept:", prob_to_accept)
            if prob_to_accept < random.random():
                #print("probability not enough")
                revert_model() # since we're using the same model we have to discard our changes
                continue
            print("accepting some randomness")

        best_prior_score = i_score
    
    print "best score", best_prior_score
    
    
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
