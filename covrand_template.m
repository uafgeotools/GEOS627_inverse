%
% covrand.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Template for HW on covariance and Gaussian random fields.
%
% calls covC.m
%

clear, close all, clc

ddir = './data/';

% load data
% m_samples : M x P
%         x : M x 1
% M is the number of model parameters describing a single sample (of a covariance matrix)
% P is the number of samples
load([ddir 'covhwdata']);
whos
break
[M,P] = size(m_samples);
xmin = min(x);
xmax = max(x);
ymin = min(m_samples(:));
ymax = max(m_samples(:));
ax0 = [xmin-5 xmax+5 ymin ymax];    % axes for plotting samples

% compute grid of distances among all points
% (D is needed for calculating the covariance matrix)
D = NaN(M,M);
for ii=1:M      % index k
    for jj=1:M  % index k'
        D(ii,jj) = abs(x(ii) - x(jj));
    end
end
% shorter way to do this: [X1,X2] = meshgrid(x,x); D = abs(X1-X2);

% random pair of points
k = randi(M,1);
kp = randi(M,1);
stlab = sprintf('d(x(%i), x(%i)) = d(%.1f, %.1f) = %.1f',k,kp,x(k),x(kp),D(k,kp));

figure; hold on;
plot(x,x*0,'k.'); plot(x([k kp]),[0 0],'ro');
xlabel('x'); title(stlab); axis tight; grid on;

% example of plotting a matrix with imagesc
figure;
imagesc(D); hold on;
axis equal, axis([-1 M+1 -1 M+1]);
plot(kp,k,'wo','markersize',10,'markerfacecolor','k');
xlabel('k'' index'); ylabel('k index');
set(gca,'xtick',[0:20:120],'ytick',[0:20:120]); colorbar
title({sprintf('Distances between pairs of points (dmax = %.1f)',max(D(:))),...
    ['EXAMPLE: ' stlab]});

%--------------------------------------------------------------------------

% START PROBLEM 3 HERE




break

% START PROBLEM 4 HERE
% CODE FOR COMPUTING NORMS OF YOUR ESTIMATED SAMPLES
% THIS ASSUMES YOU HAVE VARIABLES NAMED mc_samples and Pnew

% compute mean, std, and norm for EACH SAMPLE
mean_samples  = zeros(Pnew,1);
std_samples   = zeros(Pnew,1);
norm_samples  = zeros(Pnew,1);
for ii=1:Pnew
    mc = mc_samples(:,ii);  % sample of covariance matrix
    %mean_samples(ii)  = 
    %std_samples(ii)   = 
    %norm_samples(ii)  =    % note: Matlab's 'norm' function will not work here!
end

%==========================================================================
