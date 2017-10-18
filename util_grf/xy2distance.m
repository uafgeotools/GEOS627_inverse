function [iD,ix0,iy0,ix,iy,PA,PB] = xy2distance(nx,ny) 
%XY2DISTANCE computes the distances among all input points of a regular 2D mesh (dx = dy).
%
% The actual distances are computed as
%     D = dx*iD
% The plotting mesh is computed as
%     [X,Y] = meshgrid(xmin + ix0*dx, ymin + iy0*dx)
%
% Carl Tape, 01/2010
%

% display details and figures
idisplay = 0;

% integer index vectors
% NOTE: these start with 0 for convenience in the FFT algorithm
ix0 = [0:1:nx-1]; iy0 = [0:1:ny-1];

% generate integer ix and iy vectors
[iX,iY] = meshgrid(ix0,iy0);
ix = iX(:);
iy = iY(:);

n = nx*ny;              % number of points in 2D grid
nd = 0.5*(n^2 - n);     % number of unique distances

% indexing matrices
[PA,PB] = meshgrid([1:n],[1:n]);

% this does n^2 operations, but only nd are unique
iD2 = zeros(n,n);
kk = 0;
%disp(' i j x y');
for ii=1:nx
    xref = ix0(ii);
    for jj=1:ny
        yref = iy0(jj);
        kk = kk+1;
        iD2(kk,:) = (xref-ix).^2 + (yref-iy).^2;
        %disp(sprintf('%i %i %.2f %.2f',ii,jj,x0,y0));
    end
end

% convert to distance
iD = sqrt(iD2);

if idisplay==1
    id = iD(:);
    pA = PA(:);
    pB = PB(:);
    disp('---------------------------');
    disp(sprintf('%i (x) by %i (y) = %i gridpoints',nx,ny,n));
    disp(sprintf('%i total number of distances, %i of which are unique pairs',n^2,nd));
    for ii=1:n^2
        disp(sprintf('%i-%i (%i, %i)-(%i, %i) = %i',...
            pA(ii),pB(ii),ix(pA(ii)),iy(pA(ii)),ix(pB(ii)),iy(pB(ii)),id(ii)) );
    end
    disp('---------------------------');
end

% plot figures
if idisplay==1
    ind = [1:n]';
    ax0 = [-1 nx -1 ny];
    ax1 = [0 n+1 0 n+1];
    
    % print some output
    %whos ix iy d2 pA pB
    [ind ix+1 iy+1]
    %reshape(pA,n,n),reshape(pB,n,n), reshape(id,n,n)
    PA, PB, iD
    ud = unique(id);
    disp(sprintf('%i unique nonzero entries:',length(ud)-1));
    ud(2:end)
    
    figure; hold on;
    plot(ix,iy,'.','markersize',16);
    text(ix,iy,num2str(ind));
    axis equal, axis(ax0); grid on;
    xlabel('x (unshifted and unscaled)'); ylabel('y (unshifted and unscaled)');
    title('Indexing of points in the mesh');
    %orient tall, wysiwyg;
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); imagesc(PA); axis equal, axis(ax1);
    caxis([1 n]), colorbar
    title('Point A index');
    
    subplot(nr,nc,2); imagesc(PB); axis equal, axis(ax1);
    caxis([1 n]), colorbar
    title('Point B index');
    
    orient tall, wysiwyg;
    
    figure; nr=2; nc=1;
    subplot(nr,nc,1); imagesc(iD); axis equal, axis(ax1); colorbar
    xlabel('Index of point B'); ylabel('Index of point A');
    title(sprintf('Index distance between points A and B, max(iD) = %.2f',max(iD(:))));
    
    orient tall, wysiwyg;
    
    %-----------------------------------------------------
    % the next figures are related to an assumed covariance function
    
    % compute covariance matrix
    iL = 1;
    sigma = 1;
    stit = sprintf('(nx,ny,n) = (%i,%i,%i), iL = %i, sigma = %i',nx,ny,n,iL,sigma);
    stit2 = sprintf('%i x %i block Toeplitz with %i x %i (%i) blocks, each %i x %i',...
        n,n,nx,nx,nx*nx,ny,ny);
    R = sigma^2 * exp(-iD.^2 / (2*iL^2) );
    
    subplot(nr,nc,2); imagesc(R);
    axis equal, axis(ax1); colorbar
    xlabel('Index of point A'); ylabel('Index of point B');
    title({'Gaussian covariance between points A and B',stit2});
    
    % Gaussian sample
    A = chol(R,'lower');
    g = A*randn(n,1);
    figure; hold on;
    %pcolor( reshape(ix,ny,nx), reshape(iy,ny,nx), reshape(g,ny,nx) ); shading flat;
    imagesc( reshape(g,ny,nx) ); set(gca,'ydir','normal');
    plot(ix+1,iy+1,'o','markerfacecolor','w','markeredgecolor','b');
    text(ix+1,iy+1,num2str(ind));
    axis equal, axis(ax0+1); colorbar
    caxis(3*sigma*[-1 1]);
    xlabel('x (unshifted and unscaled)'); ylabel('y (unshifted and unscaled)');
    title(sprintf(' Cm sample: %s',stit));
    %orient tall, wysiwyg;
end

%==========================================================================
% EXAMPLE

if 0==1
    clear, clc, close all
    %nx = 3; ny = 5;      %  9 5 x 5 Toeplitz blocks
    nx = 5; ny = 3;      % 25 3 x 3 Toeplitz blocks
    xmin = -4; xmax = 12;
    Dx = xmax - xmin;
    dx = Dx/(nx-1);
    [iD,ix0,iy0,ix,iy,PA,PB] = xy2distance(nx,ny); 
    whos ix0 iy0 iD ix iy PA PB
    
    ymin = -1.5;
    xvec = xmin + ix0*dx;
    yvec = ymin + iy0*dx; 
    [X,Y] = meshgrid(xvec,yvec);
    D = dx * iD;
    % check values
    length(xvec), length(yvec), min(xvec), max(xvec), min(yvec), max(yvec)
end

%==========================================================================
