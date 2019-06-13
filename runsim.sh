#!/bin/bash

# cd /home/rihab/Documents/MATLAB/estimation/
matlab -nodisplay -r "warning('off','all');np_case=1;h=2;algo=1;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=1;h=2;algo=2;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=1;h=2;algo=3;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=1;h=3;algo=1;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=1;h=3;algo=2;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=1;h=3;algo=3;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=2;h=1;algo=1;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=2;h=1;algo=2;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=2;h=1;algo=3;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=2;h=2;algo=1;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=2;h=2;algo=2;run estimation.m;exit"
matlab -nodisplay -r "warning('off','all');np_case=2;h=2;algo=3;run estimation.m;exit"
