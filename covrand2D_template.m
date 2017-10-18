%
% covrand2D.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
% 
% Generating 2D Gaussian random fields with prescribed covariance,
% such as Gaussian or exponential.
% 
% calls
%   covC.m
%   xy2distance.m
%   xy2distance_row1.m
%

clear, clc, close all

% add path to FFT scripts for Gaussian random fields
%path(path,'/usr/local/GEOTOOLS_copy/matlab_util/util_grf/');
path(path,'./util_grf/');

%--------------------------------------------------------------------------
% USER INPUT: n, L (=iL*dx), sigma

% discretize the grid
% (assume that x and y are in units of km)
nx = 2^5; ny = nx/2; xmin = -10; xmax = 108; ymin = -20;
%nx = 2^7; ny = nx/2; xmin = -10; xmax = 108; ymin = -20;   % Problem 1-7,1-8

% correlation length of 2D Gaussian fields
% NOTE: actual correlation length is L = iL*dx
iL = 0.1250*nx;

% standard deviation of 2D fields
sigma = 0.05;

ichol = 1;          % (0,1) generate samples using Cholesky algorithm
ifourier = 0;       % (0,1) generate samples using FFT algorithm
idouble = 1;        % (0,1) =1 to remove periodic boundaries (ifourier=1 only)
iplot = 1;          % (0,1) generate plots
icov = 1;           % (1,2,3,4) type of covariance matrix (see stcovs)

% number of samples of the covariance matrix to generate
nsample = 1000;

% max n to allow for computing the n x n covariance matrix
NMAX = 2^11;

%--------------------------------------------------------------------------

n = nx*ny;      % dimension of covariance matrix
stcovs = {'Gaussian','Exponential','Circular','Matern'};
stcov = stcovs{icov};

% uniform increment in both directions (dy = dx)
Dx = xmax - xmin;
dx = Dx/(nx-1);

% grid for plotting
[iDrow1,ix0,iy0] = xy2distance_row1(nx,ny);
ymax = ymin + dx*(ny-1);
xvec = xmin + ix0*dx;
yvec = ymin + iy0*dx; 
[X,Y] = meshgrid(xvec,yvec);

% max distance between two gridpoints
% (hypotenuse of the rectangular grid)
xran = xmax-xmin;
yran = ymax-ymin;
dmax = sqrt( xran^2 + yran^2 )
%dmax = dx*sqrt((nx-1)^2 + (ny-1)^2)

% correlation length of 2D Gaussian fields
L = iL*dx;

stit0 = sprintf('L=%.2f, sigma=%.2f',L,sigma);
stit1 = sprintf('(nx,ny) = (%i,%i), L = %.2f, \\sigma = %.2f',nx,ny,L,sigma);
stit2 = sprintf('Block Toeplitz: %i x %i (%i) blocks, each %i x %i',nx,nx,nx*nx,ny,ny);

% axes limits for plotting only
axfac = 1.05;
axc = [1 n 1 n];
ax1a = [xmin xmax ymin ymax];
ax1b = [1 nx 1 ny];
axcex = axes_expand(axc,1.05,0);
ax1aex = axes_expand(ax1a,1.05,0);
ax1bex = axes_expand(ax1b,1.05,0);

% only compute the n x n matrix if n is reasonably small
if n <= NMAX
    % generate mesh of points
    [iD,ix0,iy0] = xy2distance(nx,ny);
    ymax = ymin + dx*(ny-1);
    whos dx ix0 iy0 iD

    % note: input is in indices (iD,iL)
    C = covC(iD,{icov,iL,sigma});

    % plot covariance function C(d) using a finer discretization
    dmaxt = 5*L;     % for plotting
    dcont = linspace(0,dmaxt,1000);
    ccont = covC(dcont,{icov,L,sigma});
    figure; axtemp = [0 dmaxt sigma^2*[-0.1 1.1] ];
    hold on;
    plot(axtemp(1:2),[0 0],'k');
    plot([L L],axtemp(3:4),'r');
    plot(axtemp(1:2),covC(L,{icov,L,sigma})*[1 1],'r');
    plot(dcont,ccont,'b');
    plot(iD*dx,C,'b.');       % actual values
    xlabel('Distance'); axis(axtemp);
    title(sprintf('%s covariance with L = %.2f',stcov,L));

    disp(sprintf('(nx,ny,n) = (%i,%i,%i), L = %.2f, sigma = %.3f',...
        nx,ny,n,L,sigma ));

    % vector associated with parameterization
    Atot = nx*ny*dx^2;               % NOTE: does not equal xran*yran
    dA = Atot/n;
    dAvec = ones(n,1)*dA;
    Avec = sqrt( dAvec/Atot );
    %n*dA - Atot, sum(Avec.^2)         % check

    % inverse covariance matrix for norm operations
    Cmodinv  = diag(Avec) * inv(C) * diag(Avec);            % unstable for high cond(C)
    Cdiaginv = diag(Avec) * diag(1./diag(C)) * diag(Avec);  % approximation

    %-------------------------------------------------
    % initial plots

    figure; nr=2; nc=1;
    subplot(nr,nc,1); plot(X,Y,'b.');
    axis equal, axis(ax1aex);
    xlabel('x'); ylabel('y');
    title(sprintf('Plotting grid (nx, ny, dx) = (%i, %i, %.2f), (x0, y0) = (%.1f, %.1f)',nx,ny,dx,xmin,ymin));
    subplot(nr,nc,2); plot((X-xmin)/dx,(Y-ymin)/dx,'r.');
    axis equal, axis(ax1bex-1);
    xlabel('x'); ylabel('y');
    title(sprintf('Computational grid (nx, ny) = (%i, %i)',nx,ny));
    orient tall, wysiwyg, fontsize(10)

    % plot (1) theoretical mean and (2) theoretical covariance matrix
    figure; nr=2; nc=1;
    subplot(nr,nc,1); imagesc(reshape(zeros(n,1),ny,nx) );
    set(gca,'ydir','normal');
    xlabel('x (unshifted and unscaled)'); ylabel('y (unshifted and unscaled)');
    title('theoretical mean');
    caxis(sigma*[-1 1]); colorbar;
    axis equal, axis(ax1bex);

    subplot(nr,nc,2);
    imagesc(C); shading flat;
    title({sprintf('C: %s covariance matrix with %i^2 = (%i*%i)^2 entries',stcov,n,nx,ny),stit2});
    axis equal, axis(axcex)
    caxis([0 sigma^2]); colorbar

    orient tall, wysiwyg, fontsize(11)
    
    % plot distance matrix that the covariance matrix is based on
    figure;
    imagesc(dx*iD); shading flat;
    title({sprintf('dx*iD: distance matrix with %i^2 = (%i*%i)^2 entries',n,nx,ny),stit2});
    axis equal, axis(axcex)
    caxis([0 dmax]); colorbar
    
else
    if ichol==1, error('decrease n for ichol=1, or set ichol=0'); end 

end  % if n > NMAX

%-------------------------------------------------

% ugly indexing to handle all options
if and(ichol==1, ifourier==1), izmin=1; izmax=2; end
if and(ichol==1, ifourier==0), izmin=1; izmax=1; end
if and(ichol==0, ifourier==1), izmin=2; izmax=2; end
if and(ichol==0, ifourier==0), izmin=0; izmax=0; iplot = 0; end

for iz = izmin:izmax
    
    % compute samples using Cholesky decomposition
    if and(iz==1,ichol == 1)
        disp('CHOLESKY ALGORITHM');
 
        % generate samples of the prior
        R = chol(C,'lower');
        covm_samples = zeros(n,nsample);
        for ii=1:nsample
            covm_samples(:,ii) = R * randn(n,1);
        end
    end
    
    % compute samples using a FFT method (PDF notes by Miranda Holmes)
    if and(iz==2,ifourier==1)
        
        % scale by factor of 2 to eliminate the periodic boundaries
        if idouble==1
            nxd = 2*nx; nyd = 2*ny;
        else
            nxd = nx; nyd = ny;
        end
        ix0d = [0:nxd-1];
        iy0d = [0:nyd-1];
        
        disp('FFT ALGORITHM');
        % discretization in wavenumber space
        xs = ix0d*dx;
        ys = iy0d*dx;
        kvec = k_of_x(xs);
        lvec = k_of_x(ys);
        [Kmat,Lmat] = meshgrid(kvec,lvec);
        
        % covariance function is the first row (or column)
        % note that nxd and nyd may not equal nx and ny
        %c = C(1,:);
        iDrow1 = xy2distance_row1(nxd,nyd);
        c = covC(iDrow1,{icov,iL,sigma});

        % modified version used for FFT
        cmod = reshape(c,nyd,nxd);
        cmod(nyd/2+2:end,:) = flipud(cmod(2:nyd/2,:));
        cmod(:,nxd/2+2:end) = fliplr(cmod(:,2:nxd/2));
        %figure; imagesc(cmod); axis equal; axis tight;

        % general case (using the first row of the covariance matrix)
        [kvec, lvec, CKL0] = mhfft2(xs,ys,cmod);
        CKL = real(CKL0); 
        % check that Cfunk has nominal imaginary component
        %norm(real(CKL0)), norm(imag(CKL0))
        
        % compute samples
        covm_samples = zeros(n,nsample);
        phi_h = grf2(kvec,lvec,2*CKL,nsample);  % use 2*Cfun since will take real part later
        [xss,yss,phi] = mhifft2(kvec,lvec,phi_h);
        for ii=1:nsample
            ptemp = real(phi(1:ny,1:nx,ii));   % take only the nx x ny patch
            covm_samples(:,ii) = ptemp(:);  
        end
        
    end
    
    %======================================================================
    % START YOUR PLOTS AND CALCULATIONS HERE
    % (You will also need to change the USER INPUT at the top of the code.)
    
    % plot samples
    figure; nc=2; nr=4;
    for ii = 1:min(nr*nc,nsample)
        subplot(nr,nc,ii);
        pcolor(X, Y, reshape(covm_samples(:,ii),ny,nx) ); shading flat
        axis equal, axis(ax1aex);
        caxis(3*sigma*[-1 1]);
        title(sprintf('%s sample %i/%i',stcov,ii,nsample));
    end
    xlabel('x'); ylabel('y');
    title({sprintf('%s sample %i/%i',stcov,ii,nsample),stit1});
    colorbar;
    orient tall, wysiwyg, fontsize(10)
    
    % PROBLEM 1-2
    % compute and plot the sample mean
    
   
    
    
    if n <= NMAX
        % PROBLEM 1-3
        % compute and plot covariance matrix from the samples
        
        
    end
    
    if and(iz==2,ifourier==1)
        % PROBLEM 1-7
        % vary icov and see how the GRFs change (ifourier = 1)
        % (code extracted from above)
        icvec = [1 2 3];
        nrun = length(icvec);
        nsamp = 1;
        figure;
        for xx=1:nrun
            icov = icvec(xx);
            c = covC(iDrow1,{icov,iL,sigma});
            cmod = reshape(c,nyd,nxd);
            cmod(nyd/2+2:end,:) = flipud(cmod(2:nyd/2,:));
            cmod(:,nxd/2+2:end) = fliplr(cmod(:,2:nxd/2));
            [kvec,lvec,CKL0] = mhfft2(xs,ys,cmod);
            CKL = real(CKL0); 
            
            % REPLACE THE FOLLOWING LINE WITH OTHER LINES (see grf2.m)
            phi_h = grf2(kvec,lvec,2*CKL,nsamp);
            
            [xss,yss,phi] = mhifft2(kvec,lvec,phi_h);
            ptemp = real(phi(1:ny,1:nx));
            csample = ptemp(:);  

            subplot(3,1,xx);
            pcolor(X,Y,reshape(csample,ny,nx)); shading flat; colorbar;
            axis equal, axis(ax1aex); caxis(3*sigma*[-1 1]);
            title(sprintf('%s covariance, L = %.2f, \\sigma = %.2f',stcovs{icvec(xx)},L,sigma));
        end

        % PROBLEM 1-8
        % vary L and see how the GRFs change
        

    end
    
end   % for iz

%==========================================================================
