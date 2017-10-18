function x = collocate(xmin,xmax,n)
%COLLOCATE simple collocation discretization
%
% Example: x = collocate(0,1,20)

% Aster Eqs. 1.34, 1.35
dx = (xmax-xmin)/n;
x = xmin + dx/2 + ([1:n]'-1)*dx;