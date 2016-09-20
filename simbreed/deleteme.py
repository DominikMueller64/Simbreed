#!/usr/bin/python

from time import sleep



@profile
def its_time_for_the_calculator(foo):
    """ It's time for the calculator. """
    if not isinstance(foo, int):
        return None

    a = []
    for i in range(foo):
        a.append(i)
    
    
    def so_slow(bar):
        """ Simulate a slow function. """
        sleep(5)
        return bar
    
    

    b = so_slow(a)

    c = 0
    for i in range(foo):
        c += i

    return None
