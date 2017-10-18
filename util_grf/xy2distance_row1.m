function [iDrow1,ix0,iy0] = xy2distance_row1(nx,ny) 
%XY2DISTANCE_ROW1 provides the first row of the distances among gridpoints for a 2D grid
% 
% This information will completely characterize the distances needed to
% generate stationay covariance matrices.
%
% See xy2distance.m and examples below.
%

% integer index vectors
% NOTE: these start with 0 for convenience in the FFT algorithm
ix0 = [0:1:nx-1];
iy0 = [0:1:ny-1];

% generate integer ix and iy vectors
[iX,iY] = meshgrid(ix0,iy0);
ix = iX(:)';
iy = iY(:)';

n = nx*ny;
xref = ix(1);
yref = iy(1);
iD2row1 = (xref-ix).^2 + (yref-iy).^2;
iDrow1 = sqrt(iD2row1);

%--------------------------------------------------------------------------
% EXAMPLE

if 0==1
    clear, clc, close all
    %nx = 3; ny = 5;    %  9 4 x 4 Toeplitz blocks
    nx = 5; ny = 3;     % 16 3 x 3 Toeplitz blocks
    xmin = -4; xmax = 12;
    [dx,ix0,iy0,iD,ix,iy,PA,PB] = xy2distance(xmin,xmax,nx,ny); 
    PA(1,:), PB(1,:), iD(1,:)
    
    [iDrow1,ix0,iy0] = xy2distance_row1(nx,ny);
    iDrow1
    whos iDrow1 ix0 iy0
end

%==========================================================================
