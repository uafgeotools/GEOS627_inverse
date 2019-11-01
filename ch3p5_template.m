%
% ch3p5.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Estimation of density variations using singular value decomposition.
%

clear, close all, clc

% Exercise 3.5: TSVD problem

% Load the data (d)
load('/usr/local/matlab_toolboxes/aster/cd_5.2/Exercises/chap3/prob5/ifk.mat');
whos
ndata = length(d);      % number of data
nparm = ndata;          % number of model parameters

% construct discretization vector for model (density)
ximin = 0;
ximax = 1;
xivec = collocate(ximin,ximax,nparm);

% construct discretization vector for data (gravity)
xmin = 0;
xmax = 1;
xvec = collocate(xmin,xmax,ndata);

whos xvec xivec d

% example of using plotconst_mod.m
figure; hold on;
plotconst_mod(rand(ndata,1),ximin,ximax,{'k','linewidth',2});
plotconst_mod(rand(ndata,1),ximin,ximax,{'r--','linewidth',2});

% construct design matrix G

