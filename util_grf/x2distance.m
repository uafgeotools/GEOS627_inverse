function [iD,dx,ix] = x2distance(xmin,xmax,nx) 
%X2DISTANCE compute the distances among all input points (1D)
%
% Actual grid is x = xmin + ix*dx
% Actual distances are D = dx*iD (dx > 0)
%
% Carl Tape, 01/2010
%

% generate x index vector (starts with 0)
ix = 0:1:nx-1;

% dx: actual distance increment (for plotting)
Dx = xmax - xmin;
dx = Dx/(nx-1);

% distance among all points
[X1,X2] = meshgrid(ix,ix);
iD = abs(X1-X2);
%iD2 = (X1 - X2).^2;

% plot figure
if 0==1
    nd = 0.5*(nx^2 - nx);     % number of unique distances
    n = nx;
    %figure; imagesc(iD2); axis equal, axis([0 n+1 0 n+1]); colorbar
    %xlabel('Index of point A'); ylabel('Index of point B');
    %title('Distance-squared between indices A and B');
    disp(sprintf('xmin = %.2f, xmax = %.2f, nx = %i, nd = %i',xmin,xmax,nx,nd));
    figure; imagesc(iD); axis equal, axis([0 n+1 0 n+1]); colorbar
    xlabel('Index of point A'); ylabel('Index of point B');
    title('Index distance between A and B');
end

%--------------------------------------------------------------------------
% EXAMPLE

if 0==1
    clc, close all
    nx = 16; xmin = -4; xmax = 37;
    [iD,dx,ix] = x2distance(xmin,xmax,nx);  
    
    % compute actual distances
    D = dx*iD;
    max(D(:))
    ix
    x = xmin + ix*dx
end

%==========================================================================
