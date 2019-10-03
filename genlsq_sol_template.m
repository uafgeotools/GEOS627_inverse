%
% genlsq_sol.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Demonstration of the iterative quasi-Newton method for a 4-parameter
% inversion for epicenter, origin time, and velocity.
%
% This is a simplified version of optimization.m
%
% calls forward_epicenter.m, plot_epicenters.m, plot_covsamples.m
%

clc
clear
close all
format compact
format short

%=========================================
% USER INPUT

nsamples = 1000;
irandom_initial_model = 0;      % 0(fixed), 1(random)
irandom_target_model = 0;       % 0(fixed), 1(random)
idata_errors = 2;               % 0(none),  1(random), 2(fixed)
ifig = 1;                       % 0,1

% which forward problem to test
iforward = 1;  % 1(default), 2(Tarantola 7-1)

% to print figures to PS files
bprint = false;
pdir = pwd;

% add path to Aster library (chi2inv)
addpath('/usr/local/matlab_toolboxes/aster/cd_5.2/Lib/');

%=========================================

inormalization = 1;
stnsamples = [num2str(nsamples) ' samples'];
stlabS = {'Sd(m^k)','Sm(m^k)','S(m^k) = Sd + Sm'};
      
%---------------------------------------------
% FORWARD PROBLEM    

if iforward==1
    forward_epicenter;
else
    forward_epicenter_crescent;
end

% predictions for prior and initial models (not necessary)
dprior   = d(mprior);
dinitial = d(minitial);

if ifig==1
    % plot different histograms of properties of the prior model covariance samples
    figure; nr=2; nc=2;
    for kk=1:nparm
        sigma = sigma_prior(kk);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = cov_samples_m(kk,:);
        subplot(nr,nc,kk); plot_histo(etemp,edges); ylim([0 0.4]); grid on;
        title({['mprior samples: Model parameter ' num2str(kk) ' (' mlabs{kk} ')'],
            ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ]});
    end
    if bprint, print(gcf,'-depsc',sprintf('%sgenlsq_mprior_f%i',pdir,iforward)); end

    % plot different histograms of properties of the data covariance samples
    figure; nr=4; nc=3;
    for ii=1:ndata
        sigma = sigma_obs(ii);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = cov_samples_d(ii,:);
        subplot(nr,nc,ii); plot_histo(etemp,edges); ylim([0 0.4]);
        title({['Data index ' num2str(ii)],
            ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ]});
    end
    orient tall, wysiwyg
    if bprint, print(gcf,'-depsc',sprintf('%sgenlsq_CD1_f%i',pdir,iforward)); end
    
    figure; hold on;
    plot(dobs_samples,'.-');
    p1 = plot(dprior,'bo-','linewidth',2,'markersize',10,'markerfacecolor','b','markeredgecolor','w');
    p2 = plot(dinitial,'ko-','linewidth',2,'markersize',10,'markerfacecolor','k','markeredgecolor','w');
    p3 = plot(dtarget,'ro--','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
    p4 = plot(dobs,'ro-','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
    legend([p1 p2 p3 p4],'g(mprior)','g(minitial)','g(mtarget)','g(mtarget) + errors',...
        'location','northwest');
    %title(' BLACK = d(mprior);  RED DASHED = d(mtarget);  RED = d(mtarget) + errors');
    xlim([0.5 ndata+0.5]); set(gca,'xtick',[1:ndata]);
    xlabel('Data index'); ylabel('Prediction value, g(m)'); grid on;
    if bprint, print(gcf,'-depsc',sprintf('%sgenlsq_CD2b_f%i',pdir,iforward)); end
end

%---------------------------------------------
% MISFIT FUNCTION : least squares, Tarantola (2005), Eq. 6.251
% (This calls the function d to compute the predictions.)

% data misfit
Sd = @(m,dobs,icobs) ( 0.5* ...
      transpose(d(m)-dobs) * icobs * (d(m)-dobs) );
% model misfit (related to regularization)
Sm = @(m,mprior,icprior) ( 0.5* ...
      transpose(m - mprior) * icprior * (m - mprior) );
% total misfit
S = @(m,dobs,mprior,icobs,icprior) ( Sd(m,dobs,icobs) + Sm(m,mprior,icprior) );

% initial model
%mnew = mprior;     % prior model
mnew = minitial;
Sd_0 = Sd(mnew,dobs,icobs);
Sm_0 = Sm(mnew,mprior,icprior);
S_0  = S(mnew,dobs,mprior,icobs,icprior);
stS0 = sprintf(' S(m0) = %.3f = %.3f(D) + %.3f(M)',S_0,Sd_0,Sm_0);
disp(stS0);

%///////////////////////////////
% grid search to plot misfit function
% (In many problems, the reason for using an iterativem gradient-based method 
% is that we do not have the computational resources needed to evaluate the
% misfit function over the entire model space.)

if iforward==2
    % regular grid of epicenters
    numx = 50;
    [xvec,yvec,numy,a,b,dx] = gridvec(xmin,xmax,numx,ymin,ymax);
    ng = length(xvec);
    sd0 = NaN(ng,1);
    sm0 = NaN(ng,1);
    for xx=1:length(xvec)
        mg = [xvec(xx) yvec(xx)]';
        sd0(xx) = Sd(mg,dobs,icobs);
        sm0(xx) = Sm(mg,mprior,icprior);
    end
    s0 = sd0 + sm0;
    
    % plot misfit function and its two parts: S = Sm + Sd
    X = reshape(xvec,a,b);
    Y = reshape(yvec,a,b);
    figure; nr=1; nc=3; cmin = -1e-5;
    subplot(nr,nc,1); pcolor(X,Y,reshape(s0,a,b)); shading flat;
    caxis([cmin max(s0)]); colorbar; axis equal, axis(axepi); title('S(m)');
    subplot(nr,nc,2); pcolor(X,Y,reshape(sd0,a,b)); shading flat;
    caxis([cmin max(s0)]); colorbar; axis equal, axis(axepi); title('Sd(m)');
    subplot(nr,nc,3); pcolor(X,Y,reshape(sm0,a,b)); shading flat;
    caxis([cmin max(sm0)]); colorbar; axis equal, axis(axepi); title('Sm(m)');
    % for this file: ps2eps -f misfit_f2.ps
    if bprint, print(gcf,'-dpsc',sprintf('%sgenlsq_misfit_f%i',pdir,iforward)); end
end

%///////////////////////////////

niter = input(' Select the number of iterations (< 10) or 0 to exit: ');
if niter==0, break; end

% initialize arrays
iter_vec = [0:niter]';
Sd_vec = zeros(niter+1,1);
Sm_vec = zeros(niter+1,1);
S_vec = zeros(niter+1,1);

% misfit for initial model
Sd_vec(1) = Sd_0;
Sm_vec(1) = Sm_0;
S_vec(1) = S_0;

% preconditional F (not Fhat!)
F0 = eye(nparm);

imethod = 2;    % quasi-Newton method
%imethod = 1;    % Newton method (FAILS!)
stlabels2 = {'none','newton','quasi','steepest','cg','cgpoly','vmmatrix','vmvector','srvmmatrix','srvmvector'}';
stmethod = stlabels2{imethod+1};

for nn = 1:niter
    %///////////////////////////////
    disp([' iteration ' num2str(nn) ' out of ' num2str(niter) ]);
    m     = mnew;
    delta = d(m);
    Ga    = G(m);

    % update the model: Tarantola (2005), Eq 6.319
    % (the line-search parameter is assumed to be nu = 1)
    Hhat  = icprior + Ga'*icobs*Ga;                           % approximate Hessian
    ghat  = Ga'*icobs*(delta - dobs) + icprior*(m - mprior);  % gradient
    dm    = -Hhat\ghat;     % model update
    mnew  = m + dm;         % new model

    % misfit function for new model
    % note: book-keeping only -- not used within the algorithm above
    Sd_vec(nn+1) = Sd(mnew,dobs,icobs);
    Sm_vec(nn+1) = Sm(mnew,mprior,icprior);
    S_vec(nn+1)  = S(mnew,dobs,mprior,icobs,icprior);

    disp(sprintf('%i/%i : prior, current, target:',nn,niter));
    disp([mprior mnew mtarget]);
    %///////////////////////////////
end

% misfit function values
disp('summary of misfit function:');
disp(sprintf('%8s%16s%16s%16s','iter','Sd','Sm','S = Sm + Sd'));
for nn = 1:niter+1
    disp(sprintf('%8i%16.10f%16.10f%16.10f',iter_vec(nn),Sd_vec(nn),Sm_vec(nn),S_vec(nn)));
end

if ifig==1
    % plot convergence curve
    ya = floor(min(log10([Sm_vec ; Sd_vec ; S_vec])));
    yb = ceil(max(log10([Sm_vec ; Sd_vec ; S_vec])));
    ylims = 10.^[ya yb];
    stit = [num2str(niter) ' iterations'];
    figure; hold on;
    plot(iter_vec, log10(Sd_vec),'r.-',iter_vec, log10(Sm_vec),'b.-',iter_vec, log10(S_vec),'k.-',...
        'linewidth',2,'markersize',20);
    legend(stlabS); xlim([-0.5 niter+0.5]); ylim(log10(ylims));
    set(gca,'xtick',[-1:niter+1]); grid on;
    xlabel('k, iteration'); ylabel(' log10[ S(m^k) ], misfit function'); title(stit);
    if bprint, print(gcf,'-depsc',sprintf('%sgenlsq_converge_%s_f%i',pdir,stmethod,iforward)); end
end

%///////////////////////////////
% variables needed: mpost, dpost, Gpost, cpost0, sigma_post, rho_post

% posterior model
mpost = mnew;
dpost = d(mpost);

% posterior model covariance matrix (e.g., Tarantola Eq. 3.53)
% note: cpost0 does not include normalization factors Cdfac and Cmfac
Gpost = G(mpost);
cpost0 = inv(Gpost'*icobs0*Gpost + icprior0);       % inverse of Hessian at mpost

% using matlab functions or not
if 0==1
    % posterior model uncertainties
    for i = 1:nparm, sigma_post(i) = cpost0(i,i)^(1/2); end

    % a posteriori model correlations
    rho_post = zeros(nparm,nparm);
    for i = 1:nparm
        for j = 1:nparm
            rho_post(i,j) = (cpost0(i,j)/(sigma_post(i)*sigma_post(j)));
        end
    end
else
    sigma_post = sqrt( diag(cpost0) );
    rho_post = corrcov(cpost0);
end

%///////////////////////////////

% a priori model correlations (for comparison)
rho_prior = corrcov(cprior0);

mpost, cpost0, rho_post

% posterior data covariance matrix (e.g., Tarantola Eq. 3.44)
cpost0_d = Gpost*cpost0*Gpost';
cpost0_d = (cpost0_d + cpost0_d')/2;    % force to be symmetric
sigma_post_d = sqrt( diag(cpost0_d) );  % IGNORING OFF-DIAGONAL ELEMENTS
rho_post_d = corrcov(cpost0_d);         % posterior correlation matrix
rho_prior_d = corrcov(cobs0);           % prior, for comparison

% probably cannot get a cholesky decomposition for Cpost_d
% Lpost_d = chol(cpost0_d,'lower');
% dpost_samples = zeros(ndata,nsamples);
% dcov_samples = zeros(ndata,nsamples);
% for xx=1:nsamples, randn_vecs_d(:,xx) = randn(ndata,1); end
% dcov_samples  = Lpost_d * randn_vecs_d;
% dpost_samples = repmat(dpost,1,nsamples) + dcov_samples;
% std_samples_d = std(dpost_samples');

%format long
disp(sprintf('model summary (%i iterations):',niter));
disp('    prior    initial  posterior   target');
disp([mprior minitial mpost mtarget]);
disp(sprintf('data summary (%i observations):',ndata));
disp('    prior    initial   posterior   target   actual');
disp([dprior dinitial dpost dtarget dobs]);

% Cholesky decomposition to obtain the square-root of cpost0
% NOTE: for large problems, this is not possible due to poor
%       conditioning of cpost0 or the inability to compute cpost0
Lpost = chol(cpost0,'lower');

% samples of the posterior distribution
mpost_samples = zeros(nparm,nsamples);
mcov_samples = zeros(nparm,nsamples);
for xx=1:nsamples, randn_vecs_m(:,xx) = randn(nparm,1); end
mcov_samples  = Lpost * randn_vecs_m;
mpost_samples = repmat(mpost,1,nsamples) + mcov_samples;

% compare the standard deviation with sigma_post
std_samples = std(mpost_samples');

% compare posterior model distribution with prior
% note: format statement allows for vectors (like sigma_prior)
disp('  ');
disp(' Compare model uncertainties : ');
disp(['             model parameter : ' sprintf('%13s', mlabs{:}) ]);
disp(['                       units : ' sprintf('%13s', ulabs{:}) ]);
disp(['                 sigma_prior = ' sprintf('%13.5f', sigma_prior) ]);
disp(['                  sigma_post = ' sprintf('%13.5f', sigma_post) ]);
disp(['   std(' sprintf('%6.0f', nsamples) ' mpost_samples) = ' sprintf('%13.5f', std_samples) ]);
disp(['    sigma_post / sigma_prior = ' sprintf('%13.5f', sigma_post./sigma_prior) ]);
disp('  ');

% compute the predictions associated with the posterior samples,
% then compare std_d_samples with sigma_post_d
d_samples = zeros(ndata,nsamples);
for xx=1:nsamples
    ms = mpost_samples(:,xx);
    d_samples(:,xx) = d(ms);
end
covd_samples = cov(d_samples');
rhod_samples = corrcov(covd_samples);
std_d_samples = sqrt(diag(covd_samples));
%std_d_samples = std(d_samples'); 

disp('  ');
disp(' Compare data uncertainties : ');
disp(sprintf('%16s %10s %10s %10s','prior','post','samples','post/prior'));
call = [sigma_obs(:) sigma_post_d(:) std_d_samples(:) sigma_post_d(:)./sigma_obs(:)];
for ii=1:ndata
   disp(sprintf('%6i%10.4f %10.4f %10.4f %10.4f',ii,call(ii,:)));
end
disp('  ');

if ifig==1
    % display distributions for each model parameter (nparm ROWS of cov_samples_m)
    figure; nr=2; nc=2;
    for kk=1:nparm
        sigma = sigma_post(kk);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = mcov_samples(kk,:);
        subplot(nr,nc,kk); plot_histo(etemp,edges); ylim([0 0.4]); grid on;
        stl1 = 'mpost samples';
        stl2 = ['Model parameter ' num2str(kk) ' (' mlabs{kk} ')'];
        stl3 = ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ];
        if kk==1, title({stl1,stl2,stl3});
        else title({stl2,stl3}); end
    end
    if bprint, print(gcf,'-depsc',sprintf('%sgenlsq_mpost2_hist_f%i',pdir,iforward)); end
    
    % correlation matrices and scatterplots
    %plot_covsamples(mprior_samples,rho_prior,'mprior',[],[],[],mlabs);
    plot_covsamples(mpost_samples,rho_post,'mpost',[],[],[],mlabs);
    plot_covsamples(mprior_samples,rho_prior,'mprior',mpost_samples,rho_post,'mpost',mlabs);
    if bprint
        h = length(findobj('type','figure'));
        figure(h-1); print(gcf,'-depsc',sprintf('%sgenlsq_mpost2_matrix_f%i',pdir,iforward));
        figure(h); orient tall; print(gcf,'-depsc',sprintf('%sgenlsq_mpost2_scatter_f%i',pdir,iforward));
    end
    
    % 'physical view' of estimated posterior data uncertainties
    % note: plot either sigma_post_d (from Cpost_d) or std_d_samples (from d(Cpost_samples))
    figure; hold on;
    plot(mpost_samples(1,:),mpost_samples(2,:),'c.');
    plot(mpost(1),mpost(2),'o','markersize',10,'markerfacecolor','c','markeredgecolor','w');
    %plot(mprior(1),mprior(2),'o','markersize',10,'markerfacecolor','b','markeredgecolor','w');
    %scatter(xrec,yrec,16^2,sigma_post_d,'filled','V'); title('estimated uncertainties for posterior predictions'); 
    scatter(xrec,yrec,16^2,std_d_samples,'filled','V'); title('uncertainties for posterior predictions, computed from samples'); 
    scatter(xrec,yrec,16^2,'kV');
    colormap('hot'); colorbar;
    axis equal; axis(axepi);
    %set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
    xlabel('X distance (km)'); ylabel('Y distance (km)');

    % plot solution
    % (if iforward=2, this will show a zoomed-in version as well)
    for xx=1:iforward
        % plot predictions for samples of the posterior
        plot_covsamples(d_samples,rho_post_d,'dpost',[],[],[],[]);
        
        % note: opts is set in forward_epicenter
        plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts,mpost);
        % plot the cpost0 samples and re-plot the two markers
        plot(mpost_samples(1,:),mpost_samples(2,:),'c.');
        plot(mpost(1),mpost(2),'o','markersize',10,'markerfacecolor','c','markeredgecolor','w');
        plot(mtarget(1),mtarget(2),'o','markersize',10,'markerfacecolor','r','markeredgecolor','w');
        
        title('samples of prior (blue) and posterior (cyan)');
        if xx==2, axis([-20 40 -20 50]); else axis(axepi); end
        if bprint
            print(gcf,'-depsc',sprintf('%sgenlsq_mpost_%s_f%i_v%i',pdir,stmethod,iforward,xx)); 
        end
    end
end

%==========================================================================
