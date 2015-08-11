Matlab project (including parts developed in C, MEX) to construct bisimilar symbolic abstractions of control systems and synthesize controllers employing the generated abstractions. The abstractions and controllers are stored in a BDD format, using the CUDD library.

## How to Install ##

In order to compile and install Pessoa you need execute the following commands:

`./configure --with-matlab-path=${MATLAB_ROOT_DIR}`

`make && make install`

To use pessoa add to your MATLAB path the folder /Pessoa\_1.4.
For a tutorial on how to use Pessoa visit our web: https://sites.google.com/a/cyphylab.ee.ucla.edu/pessoa/documentation/tutorial