#!/bin/bash

cd /home/rihab/Documents/MATLAB/estimation/
matlab -nodisplay -r "warning('off','all');h=1;run estimation.m;save('results1_sigma0-03_np144.mat');exit"
matlab -nodisplay -r "warning('off','all');h=2;run estimation.m;save('results1_sigma0-04_np144.mat');exit"
matlab -nodisplay -r "warning('off','all');h=3;run estimation.m;save('results1_sigma0-05_np144.mat');exit"
