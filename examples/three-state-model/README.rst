**********************************
pytram example - three state model
**********************************

This example addresses a model with three discrete states (0,1,2), i.e., two metastable
states (0,2) and a transition state (1).

The states have energy levels as given in the
file ``three-states-f_i.dat`` and, thus, stationary weights as given in the file
``three-states-pi_i.dat``. The transition probabilities are given in the file
``three-states-P_0_ij.dat``.

Further, we apply a bias (energy shifts) to increase the number of transitions when sampling
a Markov chain. The biased transition matrix is given in the file ``three-states-P_1_ij.dat``,
the full set of bias energies is given in the file ``three-states-b_K_i.dat``.

To test your pytram installation with this example, use one of the runscripts on the trajecories,
i.e., the files ``three-states-traj-0.dat`` and ``three-states-traj-1.dat``, or write your own
analysis script using the pytram API functions. Note that dTRAM requires the ``b_K_i`` values.

The estimated stationary distribution and free ebergy profile can be compared with the exact
reference values.
