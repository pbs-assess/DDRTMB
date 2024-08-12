iscam -delaydiff -mcmc 10000000 -mcsave 5000 -mcscale 10000000 -maxfn 2000 -crit 0.0001 1> runoutput.log 2>&1
iscam -delaydiff -mceval 1>> runoutput.log 2>>&1
