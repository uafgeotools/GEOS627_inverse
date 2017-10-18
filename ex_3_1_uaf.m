%
% Example 3.1
% from Parameter Estimation and Inverse Problems, 2nd edition, 2013
% by R. Aster, B. Borchers, C. Thurber
%

close all
clc
clear

cmap = 'gray';
%cmap = 'jet';
cmax = 1;

% Construct system matrix for the ray path models
% (t is used to avoid writing sqrt(2) everywhere)
t = sqrt(2);
G = [1 0 0 1 0 0 1 0 0
     0 1 0 0 1 0 0 1 0
     0 0 1 0 0 1 0 0 1
     1 1 1 0 0 0 0 0 0
     0 0 0 1 1 1 0 0 0
     0 0 0 0 0 0 1 1 1
     t 0 0 0 t 0 0 0 t
     0 0 0 0 0 0 0 0 t  ];

% Get the singular values for the system matrix
[U,S,V] = svd(G);

% Display singular values
disp('Singular values:')
diag(S)

% Find dimensions of G
[m,n] = size(G);

% Find and display system rank
disp('System rank:')
p = rank(G)

% split V into two parts: Vp and V0
Vp = V(:,1:p)

% Display null space vectors 
disp('Model null space vectors')
V0 = V(:,p+1:n)

% Display null space vectors reshaped to match tomography example geometry
disp('Model null space vectors reshaped into matrices')
m01 = reshape(V0(:,1),3,3)'
m02 = reshape(V0(:,2),3,3)'

% Display image of null space model V.,8
figure(1)
colormap(cmap)
imagesc(m01); axis equal, axis tight
colorbar, caxis([-0.5 0.5]);
set(gca,'xtick',1:3,'ytick',1:3);
xlabel('j'); ylabel('i');
display('Displaying image of null space model V.,8 (fig. 1)')
%print -deps2 c3fv8null.eps

% Display image of null space model V.,9
figure(2)
colormap(cmap)
imagesc(m02); axis equal, axis tight
colorbar, caxis([-0.5 0.5]);
set(gca,'xtick',1:3,'ytick',1:3);
xlabel('j'); ylabel('i');
display('Displaying image of null space model V.,9 (fig. 2)')
%print -deps2 c3fv9null.eps

% Display data space null vector
disp('Data null space vector')
U0 = U(:,p+1)

% Find and display model resolution matrix
Rm = Vp*Vp';    % note: this is NOT the identity matrix

% Display full model resolution matrix and colorbar
figure(3)
colormap(cmap)
imagesc(Rm); axis equal, axis tight
colorbar, caxis([0 cmax])
set(gca,'xtick',1:9,'ytick',1:9);
xlabel('k'''); ylabel('k');
display('Displaying full model resolution matrix (fig. 3)')
%print -deps2 c3fR.eps

Rmdiag = reshape(diag(Rm),3,3)'
% Display reshaped diagonal elements of the model resolution matrix and colorbar
figure(4)
colormap(cmap)
imagesc(Rmdiag); axis equal, axis tight
colorbar, caxis([0 cmax]);
set(gca,'xtick',1:3,'ytick',1:3);
xlabel('j'); ylabel('i');
display('Displaying diagonal elements of model resolution matrix (fig. 4)')
%print -deps2 c3fRdiag.eps

% Spike resolution test

% Construct spike model
mtest = zeros(n,1);
mtest(5) = 1;

% Get noise free data for the spike model (forward problem)
dtest = G*mtest;

% Display spike model and noise free data
disp('model spike and predicted data')
mtest
dtest

% Display recovered model from spike test
disp('recovered model from the spike test')
% Let Matlab compute the pseudoinverse for us with pinv
mdagger = pinv(G)*dtest

% Display reshaped noise free spike model with color bar
figure(5)
colormap(cmap)
imagesc(reshape(mtest,3,3))'; axis equal, axis tight
colorbar, caxis([0 cmax])
set(gca,'xtick',1:3);
set(gca,'ytick',1:3);
xlabel('j'); ylabel('i');
display('Displaying noise free spike test model (fig. 5)')
%print -deps2 c3fspike.eps

% Display reshaped recovered spike model with colorbar
figure(6)
colormap(cmap)
imagesc(reshape(mdagger,3,3)'); axis equal, axis tight
colorbar, caxis([0 cmax])
set(gca,'xtick',1:3);
set(gca,'ytick',1:3);
xlabel('j'); ylabel('i');
display('Displaying recovered model for spike test (fig. 6)')
%print -deps2 c3fspike_recov.eps
