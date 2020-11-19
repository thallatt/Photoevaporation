#! /usr/bin/env gnuplot 

set output "~/Documents/Work/MSc/code/evaporation.eps"
load("~/Documents/Work/MSc/code/plotting/settings.gnu")
set lmargin at screen 0.13
set rmargin at screen 0.95
set tmargin at screen 0.96
set bmargin at screen 0.12


set xrange [1:3e2]
set yrange [1e-2:5]

#set mytics 5
set mxtics 10

set format y "10^{%L}"
set format x "10^{%L}" 

set xlabel "Period [days]" offset -1.5,-1 font ",30"
set ylabel "Number of Planets" offset -4.5,0 font ",30"

set logscale y 10
set logscale x 10

set label 1 at gr 0.45, 0.9 "P+18, sub-Saturns" tc rgb darkred font ",30" rotate by 0 front
set arrow 1 from gr 0.85, 0.89 to gr 0.94,0.87 lw 3 lc rgb darkred dt 1 front
#set label 2 at gr 0.04, 0.33 "m^{-1.5}" tc rgb darknavy font ",19" rotate by 0 front
#set arrow 2 from gr 0.08, 0.31 to gr 0.16,0.29 lw 3 lc rgb darknavy dt 1 front
#set label 3 at gr 0.04, 0.42 "m^{-1}" tc rgb lightnavy font ",19" rotate by 0 front
#set arrow 3 from gr 0.08, 0.4 to gr 0.16,0.38 lw 3 lc rgb lightnavy dt 1 front
#set label 4 at gr 0.04, 0.51 "m^{-0.5}" tc rgb lightblue font ",19" rotate by 0 front
#set arrow 4 from gr 0.08, 0.49 to gr 0.16,0.47 lw 3 lc rgb lightblue dt 1 front
#set label 5 at gr 0.32, 0.8 "lognormal" tc rgb midorange font ",19" rotate by 0 front
#set arrow 5 from gr 0.39, 0.77 to gr 0.39,0.71 lw 3 lc rgb midorange dt 1 front
#set label 6 at gr 0.34, 0.15 "Rayleigh" tc rgb darkgreen font ",19" rotate by 0 front
#set arrow 6 from gr 0.33, 0.15 to gr 0.25,0.15 lw 3 lc rgb darkgreen dt 1 front

plot "/Users/Tim/Documents/Work/MSc/code/evaporation_rayleigh_mean=13.0_output_HD_R_8.csv" using ($2):(1*$3) with lines dt 3 lw 10 lc rgb darkblue title "",\
"/Users/Tim/Documents/Work/MSc/code/evaporation_lognormal_mean=2.7_dev=0.5_output_HD_R_8.csv" using ($2):(1*$3) with lines dt 3 lw 10 lc rgb midorange title "",\
"~/Documents/Work/MSc/petigura_saturns.csv" using ($1):($2):($4):($3) with errorbars pointtype 7 pointsize 2 lw 5 lc rgb darkred title ""
#"/Users/Tim/Documents/Work/MSc/code/evaporation_powerlaw_index=1.5_output.csv" using ($2):(1*$3) with lines dt 3 lw 5 lc rgb darknavy title "",\
#"/Users/Tim/Documents/Work/MSc/code/evaporation_powerlaw_index=1.0_output.csv" using ($2):(1*$3) with lines dt 3 lw 5 lc rgb lightnavy title "",\
#"/Users/Tim/Documents/Work/MSc/code/evaporation_powerlaw_index=0.5_output.csv" using ($2):(1*$3) with lines dt 3 lw 5 lc rgb lightblue title "",\
#"/Users/Tim/Documents/Work/MSc/code/evaporation_rayleigh_mean=10.0_output.csv" using ($2):(1*$3) with lines dt 3 lw 5 lc rgb lightgreen title "",\
#"/Users/Tim/Documents/Work/MSc/code/evaporation_lognormal_mean=2.995732273553991_dev=0.3_output.csv" using ($2):(1*$3) with lines dt 3 lw 5 lc rgb midorange title ""
#
