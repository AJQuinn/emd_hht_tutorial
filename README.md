# Time-Frequency analysis with Empirical Mode Decomposition

This tutorial introduces two concepts:
1) High-resolution time-frequency analysis using Instantaneous Frequency
estimation
2) Separation of a complicated signal into simpler oscillatory
components.

Taken together, these steps allow us to compute time-frequency plots at a
much higher resolution than available with most standard methods.

Here, we will analyse at some simulations and some example data from the
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

Finally, the ```simulate.m``` function provided here is an abridged version of the simulation
functions in the CFC toolbox. For more details see:

https://github.com/AJQuinn/cfc

# The Tutorial

To run the tutorial, open emd_time_frequency_tutorial.m in MatLab and run
through the script one-cell at a time.
