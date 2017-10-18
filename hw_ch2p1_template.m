%
% hw_ch2p1.m
%
% from Parameter Estimation and Inverse Problems, 2nd edition, 2013
% by R. Aster, B. Borchers, C. Thurber
%

close all, clc, clear

M = 6;  % number of data points
N = 2;  % number of model parameters

% vector of distances
% start of refraction
xstart = 5+1/3;
% end of refraction
xmax = 5*xstart;
x = linspace(ceil(xstart),xmax,M)';

% vector of first arrival times
t = [ 3.4935 4.2853 5.1374 5.8181 6.8632 8.1841 ]';


