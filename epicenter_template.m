%
% epicenter.m
% Tarantola 2005, Problem 7.1
%
% Be sure to understand lab_sampling.m before starting this problem.
%
% Carl Tape 2013-01-30
%

clear, close all, clc

% ENTER DATA AND OTHER PARAMETERS
dx = 0.2;       % grid spacing, km


% travel time between any two points
t = @(x,y,xr,yr) ( YOUR FUNCTION HERE );

% N x 1 vector of predicted traveltimes
tm = @(x,y) ([  t(x,y,xr(1),yr(1))
                t(x,y,xr(2),yr(2))
                t(x,y,xr(3),yr(3))
                t(x,y,xr(4),yr(4))
                t(x,y,xr(5),yr(5))
                t(x,y,xr(6),yr(6)) 
                   ]);
               
               
% range of model space (Tarantola Figure 7.1)
% note: length unit is km
xmin = 0; xmax = 22;
ymin = -2; ymax = 30;
xvec = xmin:dx:xmax; nx = length(xvec);
yvec = ymin:dx:ymax; ny = length(yvec);
% X: ny x nx matrix of x values (for plotting)
% Y: ny x nx matrix of y values (for plotting)
[X,Y] = meshgrid(xvec,yvec);
% x: ng x 1 vector of x values (for function evaluation)
% y: ng x 1 vector of y values (for function evaluation)
x = X(:);
y = Y(:);
ng = length(x);     % number of gridpoints
dA = dx^2;          % area of each cell
disp(sprintf('model space is %i x %i = %i',nx,ny,ng));
disp(sprintf('model cell is %.2f km x %.2f km = %.4f km^2',dx,dx,dA));

% misfit function, given a SINGLE model m=(x,y) and a 6 x 1 vector of data
F = @(x,y,tobs) ( YOUR FUNCTION HERE );


%==========================================================================
