# Automatic_RFIC_Synth
Repository of codes used to automatically configure a TX side RF circuit

***Spectre_Single_Point.py***
Using spectre_single_point.py, we can specify a set of input parameters to the circuit, invoke spectre through the CLI, run multiple analyses and print the outputs on the console

invoked using *python spectre_single_point.py*

***Ckt_Optimizer.py***
Using ckt_optimizer.py, we can specify a set of output constraints to the PA circuit and the circuit will use gradient descent algorithm to find the optimal circuit parameters satisfying the specified constraints.

invoked using *python ckt_optimizer.py | tee logfile*
