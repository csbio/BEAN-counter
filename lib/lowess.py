# Modified version of agramfort's lowess implementation, located at 
# https://gist.github.com/agramfort/850437
# and licensed under the BSD (3-clause)

# BSD-3-Clause
# Copyright 2018 Regents of the University of Minnesota
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np
import os, gzip, cPickle
from scipy.linalg import solve
from multiprocessing.pool import Pool
from contextlib import closing

def py_lowess(x, y, f=2. / 3., iter=3, num_cores = 2):
    '''
    Performs a local scatterplot smoothing on y given x (LOWESS).
    Results line up with MATLAB's 'smooth' function.
    '''

    # Sets globals to enable efficient threading
    global x_
    global y_
    global h
    global delta
    global starts
    global stops
    global r

    # First ensures that the input arrays are properly formatted
    if len(x) != len(y): return np.nan 
    if len(x) == 0 or len(y) == 0: return np.nan
    if isinstance(x, np.ndarray):
        if [i > 1 for i in x.shape].count(True) > 1: return np.nan
        else: x = np.reshape(x.flatten(), max(x.shape))
    if isinstance(y, np.ndarray):
        if [i > 1 for i in y.shape].count(True) > 1: return np.nan
        else: y = np.reshape(y.flatten(), max(y.shape)) 

    # Sorting by x enables some more efficient schemes
    sort_inds = x.argsort()
    unsort_inds = sort_inds.argsort()
    x = x[sort_inds]
    y = y[sort_inds]

    # Giving the globals their values **after** x and y are sorted
    # Otherwise bad things happen (incorrect results, singular matrices)
    x_ = x
    y_ = y

    # Calculate size of span
    n = len(x)
    r = int(np.ceil(f * n))

    # Calculate each x value's span (distance on the x axis from x)
    if num_cores > 1:
        with closing(Pool(processes = num_cores)) as pool:
            h = np.array(pool.map(calc_span, range(n)))
            pool.terminate()
    else:
        h = np.array(map(calc_span, range(n)))

    # Calculate the start and end span indices for each x value
    starts, stops = span_inds(x, h, r)

    # Calculate the y estimates, iterate the specified number of times
    delta = np.ones(n)
    for iteration in range(iter):
        #print 'iteration:', iteration
        if num_cores > 1:
            with closing(Pool(processes=num_cores)) as pool:
                yest = np.array(pool.map(calc_yest, range(n)))
                pool.terminate()
        else:
            yest = np.array(map(calc_yest, range(n)))
        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2
    return yest[unsort_inds]

def cube(x):
    '''
    This is much faster than x ** 3, using x ** 3 just happens
    to be a bottleneck!
    '''
    return x * x * x

def span_inds(x, h, r):
    '''
    Calculates span indices for each value in sorted 1D array x,
    given the span in units of distance on the x axis (h). r
    supplies the default width of the span in indices (before
    adjustment for ties, etc.).
    '''
    assert np.all(np.diff(x) >= 0), "Argument x must be sorted in ascending order!"
    n = len(x)
    starts_ = np.zeros_like(x, dtype = 'int')
    stops_ = np.zeros_like(x, dtype = 'int')
    start = 0
    stop = r
    for i,foc_x in enumerate(x):
        while abs(foc_x - x[start]) > h[i]:
            start += 1
        if stop < n:
            while abs(foc_x - x[stop]) < h[i]:
                stop += 1
                if stop == n:
                    break

        starts_[i] = start
        stops_[i] = stop - 1
    return starts_, stops_

def calc_span(i):
    global x_
    global r
    return np.partition(np.abs(x_ - x_[i]), r)[r]

def calc_yest(i):
    #x_i, x_span, y_span, h_i, delta_span = data
    # Ensure that we're referencing global variables
    global x_
    global y_
    global h
    global delta
   
    # If this has already converged, as evidenced by all delta values
    # in the span being zero, then just return the previous y estimate.
    # NOTE: this completely depends on calculating the span such that
    # the interval is exclusive. This way, it is guaranteed (hopefully)
    # that all values in w are nonzero and that the only thing that can
    # cause all weights to be zero is if all delta values within the
    # span are zero.
    delta_span = delta[starts[i]:(stops[i] + 1)]
    if np.allclose(delta_span, 0):
        return y_[i]
    
    x_span = x_[starts[i]:(stops[i] + 1)]
    y_span = y_[starts[i]:(stops[i] + 1)]

    w = np.clip(np.abs((x_span - x_[i]) / h[i]), 0.0, 1.0)
    w = cube(1 - cube(w))
    weights = delta_span * w
    b = np.array([np.sum(weights * y_span), np.sum(weights * y_span * x_span)])
    A = np.array([[np.sum(weights), np.sum(weights * x_span)],
                  [np.sum(weights * x_span), np.sum(weights * x_span * x_span)]])
    beta = solve(A, b)
    yest_i = beta[0] + beta[1] * x_[i]
    return yest_i
