# Plots the observed counting process (ocp) associated with a spike train
# Change the value of sampling_frequency to suit your data! 
sampling_frequency=15000.
set xlabel "Time (s)"
set ylabel "Cumulative spike number"
set grid
unset key
s=0
plot "-" u ($1/sampling_frequency):(s=s+1) with steps
