SUM-PRODUCT Implementation
--------------------------

Author: Mike Schachter (mike.schachter@gmail.com)


REQUIREMENTS
------------
python 2.5 or higher
numpy (most versions will do)
networkx package: (easy_install networkx==1.5)


HOW TO RUN (Linux and Mac OSX)
------------------------------
Running problem 3b is as simple as:

export PYTHONPATH=$PYTHONPATH:.
python example.py

The marginals are printed out at the end. They should be:

P[x1] = (0.59, 0.41)
P[x2] = (0.10, 0.90)
P[x3] = (0.60, 0.40)
P[x4] = (0.07, 0.93)
P[x5] = (0.57, 0.43)
P[x6] = (0.13, 0.87)

HOW TO RUN (Windows)
--------------------

Stop using windows! Go find a Linux or Mac machine.


OVERVIEW
--------

This is an implementation of SUM-PRODUCT that works using
a synchronous parallel flooding schedule. It's not completely
local or completely parallel, but could easily be made as such.

The bulk of the code is in sum_product.py. There are several
important classes:

UGPGraph: contains the graph structure and potentials

Message: a data structure for containing nested sets of messages

SumProduct: a class for running the sum-product algorithm with the
flooding schedule.

SPNode: a helper data-structure for the SumProduct Class. 

See the source code for more documentation.

