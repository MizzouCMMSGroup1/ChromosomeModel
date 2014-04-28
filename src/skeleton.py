
import numpy

# bring in our cleaned data
IF_FILENAME = "if_data_stripped.csv"
if_data = numpy.loadtxt(IF_FILENAME, delimiter=',')

# weighting scores from paper, we are using megabase resolution

W1 = 1.0
W2 = 1.5
W3 = 1.5
W4 = 1.5

minimum_threshold_squared = 7 # micrometers

def contact_score():
    '''
    minimize the distance (but keep above min_threshold) between non-sequential pairs that have affinity
    '''
    return 1


def noncontact_score():
    '''
    maximize the distance (but keep below max_threshold) between non-sequential pairs that don't have affinity
    '''
    return 1


def pair_smoothing():
    '''
    keep adjacent contacts (eg |i-j|==1) with slightly lower score than above so they are prioritized for optimization
    '''
    return 1


def model_score(model=model):
    return contact_score() + noncontact_score() + pair_smoothing()


def main():
    
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
