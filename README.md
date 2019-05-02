# Time-Frequency analysis with Empirical Mode Decomposition

This tutorial introduces two concepts:
1) High-resolution time-frequency analysis using Instantaneous Frequency
estimation
2) Separation of a complicated signal into simpler oscillatory
components.

Taken together, these steps allow us to compute time-frequency plots at a
much higher resolution than available with most standard methods.

Here we will analyse at some simulations and some example data using the
fieldtrip toolbox

# Getting started

To run this tutorial you will need to  be running matlab R2018a or more
recent with the signal processing toolbox.
The  Empirical Mode Decomposition functions used here are not available
in older versions.


Next, you will need a recent version of Fieldtrip on the  matlab path.
Fieldtrip and a host of tutorials are available here:

http://www.fieldtriptoolbox.org

Fieldtrip is used to compute a Hanning-Taper time-frequency transform for
comparison with the Empirical Mode Decomposition.

Finally, you will  need a copy of a Matlab Cross-Frequency Coupling toolbox
which can be downloaded here:

https://github.com/AJQuinn/cfc
git clone https://github.com/AJQuinn/cfc.git

We will use the simulation functions in this toolbox to generate some
signals to illustrate the EMD.
