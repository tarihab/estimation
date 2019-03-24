%% Code for estimation of scalar field

% Domain is unit square
xborder = [0 1 1 0];
yborder = [0 0 1 1];

nagents = 5;

p = 100;  % no.of parameters
xc = 0.05:0.1:1;  % x-coordinates of RBF centres
yc = 0.05:0.1:1;  % y-coordinates of RBF centres
sigma = 0.03;  % std. deviation of RBFs
