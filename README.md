# Monte-Carlo-Particle-Simulation

## Code Description
### Input
    The input file is structured as follows:
        <function name> <trials> <stick length> <iterations>
        - Function name calls for either "circle" or "buffon" estimation
        - Trials determines how many random samples are taken each iteration, accepts integers
        - Stick length is only relevant for the Buffon's needle test where it defines how long the thrown stick is, accepts floats.
        - Iterations defines how many times the test is run with the same parameters
    The input file is benchmark.txt
### Functions
    randomVal
        - Creates a random value between 0 and 1 or optionally self-selected parameters min and max that are sourced from the mt19937 generator. 
    circleSampling
        - Runs the circle estimation algorithm for a given amount of iterations. Works by calculating x^2 + y^2 and checking whether it is less than 1.
    buffonNeedle
        - Runs the Buffon's needle estimation algorithm for a given amount of iterations and a given length of stick (needle). Default length is 1. The width of planks is set to 1. 
    timing
        - This function runs the function passed to it and times its execution and returns the results in a Result struct. 
    statistics
        - This function calculates all the parameters of the function iteration including mean, standard deviation, total running time, figure-of-merit, and the JB normality estimate. It returns a Stats struct where all this information is stored. 
    printStats
        - This function prints out the information stored inside the Stats struct.
    
    The main function loop reads the benchmark.txt file and accordingly calls any functions there. Once each line of the benchmark is run, the results are printed out.

## Completed bonus exercises

1. Statistical test for normality
    - This was achieved by implementing a Jarque-Bera statistical test for defining the normality of a distribution. The 3rd and 4th moment terms are calculated in the same method as the standard deviation and the test statistic JB is assigned to the Stats struct. The test shows that once the amount of trials per iteration grows high enough, the test is passed (JB < 5.991). 

## Performance
The best prediction result was achieved (of the ones tested) with a needle length of 2, 10 000 000 trials and 10 iterations. This returned a FOM of around 400k. 
On a M2 MacBook, this function took 12.9 seconds to run. 

## Time usage
Assignment 1 : 5h

## Misc
Compiled with
 >g++ -std=c++20 -g -Wall -Wextra -o mc mc.cpp