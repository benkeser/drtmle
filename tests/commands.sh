#!/bin/bash

# lines to scp files to snail and rhino 

# get to directory with simulation files
cd ~/Dropbox/Dissertation/one\ step\ dr\ inference/drtmle/tests

# scp cent and sce
scp centBiv.R sceBiv.sh centSix.R sceSix.sh dbenkese@snail.fhcrc.org:~/dral

# log-in to snail
ssh dbenkese@snail.fhcrc.org
cd dral
scp centBiv.R sceBiv.sh centSix.R sceSix.sh dbenkese@rhino.fhcrc.org:~/dral
ssh dbenkese@rhino.fhcrc.org
cd dral
./sceBiv.sh ./centBiv.R newTwo_adaptive
./sceSix.sh ./centSix.R newSix





# get results to snail
scp ~/dral/out/newSim6_allOut.RData dbenkese@snail.fhcrc.org:~/dral
scp ~/dral/out/newSim2_allOut.RData dbenkese@snail.fhcrc.org:~/dral

# get results to computer
cd ~/Dropbox/Dissertation/one\ step\ dr\ inference/computerCode/computerOut
scp dbenkese@snail.fhcrc.org:~/dral/newSim2_allOut.RData . 
scp dbenkese@snail.fhcrc.org:~/dral/newSim6_allOut.RData . 