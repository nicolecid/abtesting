from scipy import stats
from scipy.stats import t as t_dist
from scipy.stats import chi2
import math
from abtesting_test import *

# You can comment out these lines! They are just here to help follow along to the tutorial.
#print(t_dist.cdf(-2, 20)) # should print .02963
#print(t_dist.cdf(2, 20)) # positive t-score (bad), should print .97036 (= 1 - .2963)
#
#print(chi2.cdf(23.6, 12)) # prints 0.976
#print(1 - chi2.cdf(23.6, 12)) # prints 1 - 0.976 = 0.023 (yay!)

# TODO: Fill in the following functions! Be sure to delete "pass" when you want to use/run a function!
# NOTE: You should not be using any outside libraries or functions other than the simple operators (+, **, etc)
# and the specifically mentioned functions (i.e. round, cdf functions...)

def slice_2D(list_2D, start_row, end_row, start_col, end_col):
    '''
    Splices a the 2D list via start_row:end_row and start_col:end_col
    :param list: list of list of numbers
    :param nums: start_row, end_row, start_col, end_col
    :return: the spliced 2D list (ending indices are exclsive)
    '''
    to_append = []
    for l in range(start_row, end_row):
        to_append.append(list_2D[l][start_col:end_col])

    return to_append

def get_avg(nums):
    '''
    Helper function for calculating the average of a sample.
    :param nums: list of numbers
    :return: average of list
    '''
    length = len(nums)
    adding = sum(nums)
    average = adding / length
    return average

def get_stdev(nums):
    '''
    Helper function for calculating the standard deviation of a sample.
    :param nums: list of numbers
    :return: standard deviation of list
    '''
    average = get_avg(nums)
    length = len(nums) - 1
    sumation = 0
    for i in nums: 
        substraction = i - average
        square = substraction**2
        sumation += square
    div = sumation / length
    standev = div**(1/2)
    return standev

def get_standard_error(a, b):
    '''
    Helper function for calculating the standard error, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: standard error of a and b (see studio 6 guide for this equation!)
    '''
    std_a = get_stdev(a)
    std_b = get_stdev(b)
    n_a = len(a)
    n_b = len(b)
    std_a_sqr = std_a**2
    std_b_sqr = std_b**2
    div_a = std_a_sqr/n_a
    div_b = std_b_sqr/n_b
    summa = div_a + div_b
    stand_error = summa**(1/2)
    
    return stand_error

def get_2_sample_df(a, b):
    '''
    Calculates the combined degrees of freedom between two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: integer representing the degrees of freedom between a and b (see studio 6 guide for this equation!)
    HINT: you can use Math.round() to help you round!
    '''
    std_a = get_stdev(a)
    std_b = get_stdev(b)
    n_a = len(a)
    n_b = len(b)
    std_a_sqr = std_a**2
    std_b_sqr = std_b**2
    div_a = (std_a_sqr/n_a)**2
    div_b = (std_b_sqr/n_b)**2
    standard_error = get_standard_error(a,b)
    stand_error_4 = standard_error**4
    div_n_a = div_a/(n_a - 1)
    div_n_b = div_b/(n_b - 1)
    summa = div_n_a + div_n_b
    to_round = stand_error_4 / summa
    df = round(to_round)
    return df

def get_t_score(a, b):
    '''
    Calculates the t-score, given two samples.
    :param a: list of numbers
    :param b: list of numbers
    :return: number representing the t-score given lists a and b (see studio 6 guide for this equation!)
    '''
    avg_a = get_avg(a)
    avg_b = get_avg(b)
    subs = avg_a - avg_b
    stand_error = get_standard_error(a, b)
    t_score = subs/stand_error
    if t_score > 0:
        t_score = -t_score

    return t_score

def perform_2_sample_t_test(a, b):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates a p-value by performing a 2-sample t-test, given two lists of numbers.
    :param a: list of numbers
    :param b: list of numbers
    :return: calculated p-value
    HINT: the t_dist.cdf() function might come in handy!
    '''
    t_score = get_t_score(a,b)
    df = get_2_sample_df(a,b)
    p_value = t_dist.cdf(t_score, df)
    return p_value


# [OPTIONAL] Some helper functions that might be helpful in get_expected_grid().
def row_sum(observed_grid, ele_row):
    summa = sum(observed_grid[ele_row])
    return summa
    
def col_sum(observed_grid, ele_col):
    summa = 0
    n_rows = len(observed_grid)
    for i in range(n_rows):
        summa += observed_grid[i][ele_col]
    return summa

def total_sum(observed_grid):
    tot_sum = 0
    n_rows = len(observed_grid)
    n_cols = len(observed_grid[0])
    
    for i in range(n_rows):
        for j in range(n_cols):
            tot_sum+= observed_grid[i][j]
    return tot_sum
    
def calculate_expected(row_sum, col_sum, tot_sum):
    row_col = row_sum * col_sum
    div = row_col / tot_sum
    return div

def get_expected_grid(observed_grid):
    '''
    Calculates the expected counts, given the observed counts.
    ** DO NOT modify the parameter, observed_grid. **
    :param observed_grid: 2D list of observed counts
    :return: 2D list of expected counts
    HINT: To clean up this calculation, consider filling in the optional helper functions below!
    '''
    grid = []
    n_rows = len(observed_grid)
    n_cols = len(observed_grid[0])
    tot_sum = total_sum(observed_grid)
    
    for i in range(n_rows):
        grid_row = []
        for j in range(n_cols):
            sum_row = row_sum(observed_grid, i)
            sum_col = col_sum(observed_grid, j)
            expected_val = calculate_expected(sum_row, sum_col, tot_sum)
            grid_row.append(expected_val)
        grid.append(grid_row)
    return grid
    
def df_chi2(observed_grid):
    '''
    Calculates the degrees of freedom of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: degrees of freedom of expected counts (see studio 6 guide for this equation!)
    '''
    n_rows = len(observed_grid)
    n_cols = len(observed_grid[0])
    
    df = (n_rows-1)*(n_cols-1)
    return df

def chi2_value(observed_grid):
    '''
    Calculates the chi^2 value of the expected counts.
    :param observed_grid: 2D list of observed counts
    :return: associated chi^2 value of expected counts (see studio 6 guide for this equation!)
    '''
    grid = get_expected_grid(observed_grid)
    n_cols = len(observed_grid[0])
    n_rows = len(observed_grid)
    chi = [[0 for i in range(n_cols)] for j in range(n_rows)]
    
    for i in range(n_rows):
        for j in range(n_cols):
            expected = grid[i][j]
            obs_grid = observed_grid[i][j]
            subs = (expected - obs_grid)**2
            chi[i][j] = subs/expected
    tot_sum = total_sum(chi)
    return tot_sum

def perform_chi2_homogeneity_test(observed_grid):
    '''
    ** DO NOT CHANGE THE NAME OF THIS FUNCTION!! ** (this will mess with our autograder)
    Calculates the p-value by performing a chi^2 test, given a list of observed counts
    :param observed_grid: 2D list of observed counts
    :return: calculated p-value
    HINT: the chi2.cdf() function might come in handy!
    '''
    val = chi2_value(observed_grid)
    df = df_chi2(observed_grid)
    p_value = (1- chi2.cdf(val,df))
    return p_value

# These commented out lines are for testing your main functions. 
# Please uncomment them when finished with your implementation and confirm you get the same values :)
def data_to_num_list(s):
  '''
    Takes a copy and pasted row/col from a spreadsheet and produces a usable list of nums. 
    This will be useful when you need to run your tests on your cleaned log data!
    :param str: string holding data
    :return: the spliced list of numbers
    '''
  return list(map(float, s.split()))


## t_test 1:
#a_t1_list = data_to_num_list(a1)
#b_t1_list = data_to_num_list(b1)
#print(get_t_score(a_t1_list, b_t1_list)) # this should be -129.500
#print(perform_2_sample_t_test(a_t1_list, b_t1_list)) # this should be 0.0000
## why do you think this is? Take a peek at a1 and b1 in abtesting_test.py :)
#
#
## t_test 2:
#a_t2_list = data_to_num_list(a2)
#b_t2_list = data_to_num_list(b2)
#print(get_t_score(a_t2_list, b_t2_list)) # this should be -1.48834
#print(perform_2_sample_t_test(a_t2_list, b_t2_list)) # this should be .082379
#
## t_test 3:
#a_t3_list = data_to_num_list(a3)
#b_t3_list = data_to_num_list(b3)
#print(get_t_score(a_t3_list, b_t3_list)) # this should be -2.88969
#print(perform_2_sample_t_test(a_t3_list, b_t3_list)) # this should be .005091

"""
# chi2_test 1:
a_c1_list = data_to_num_list(a_count_1) 
b_c1_list = data_to_num_list(b_count_1)
c1_observed_grid = [a_c1_list, b_c1_list]
print(chi2_value(c1_observed_grid)) # this should be 4.103536
print(perform_chi2_homogeneity_test(c1_observed_grid)) # this should be .0427939


# chi2_test 2:
a_c2_list = data_to_num_list(a_count_2) 
b_c2_list = data_to_num_list(b_count_2)
c2_observed_grid = [a_c2_list, b_c2_list]
print(chi2_value(c2_observed_grid)) # this should be 33.86444
print(perform_chi2_homogeneity_test(c2_observed_grid)) # this should be 0.0000
# Again, why do you think this is? Take a peek at a_count_2 and b_count_2 in abtesting_test.py :)

# chi2_test 3:
a_c3_list = data_to_num_list(a_count_3) 
b_c3_list = data_to_num_list(b_count_3)
c3_observed_grid = [a_c3_list, b_c3_list]
print(chi2_value(c3_observed_grid)) # this should be .3119402
print(perform_chi2_homogeneity_test(c3_observed_grid)) # this should be .57649202
"""


