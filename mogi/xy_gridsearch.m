% xy_gridsearch.m

clear, close all, clc

ddir = './';
fsize = 14;
msize = 18;

% read the interferogram, then plot it
read_data;
whos obs_rngchg
[xvec,yvec] = plot_model(obs_rngchg);

iokay = find(isnan(obs_rngchg)==0);

%==========================================================================
% FIX Z AND dV, SEARCH OVER X AND Y

% grid for 2D search
xs = [19:0.2:22];
ys = [21:0.2:23];
zs = -2.58;
V = 0.0034;

nx = length(xs);
ny = length(ys);
ng = nx*ny;

disp(sprintf('fixed zs = %.2f km, V = %.4f, searching over (xs,ys)',zs,V));
misfit = NaN(nx,ny);
kk = 0;
for k = 1:nx
    for l = 1:ny
        kk = kk+1;
        syn_rngchg = mogi2insar(xs(k),ys(l),zs,V,0,0);  
        misfit(k,l) = sum((obs_rngchg(iokay) - syn_rngchg(iokay)).^2);
    end
    disp(sprintf('source %i/%i is (xs = %.2f km, ys = %.2f km)',kk,ng,xs(k),ys(l)));
end

% get the minimum
[indx,indy] = find(misfit == min(misfit(:)));
disp(sprintf('Source X coordinate: %.2f km',xs(indx)));
disp(sprintf('Source Y coordinate: %.2f km',ys(indy)));

% plot cross section of misfit function
figure;
imagesc(xs,ys,misfit');
set(gca,'ydir','normal'); hold on;
plot(xs(indx),ys(indy),'kp','markersize',msize,'markerfacecolor','w');
colorbar
set(gca,'fontsize',fsize-2);
xlabel('Easting [km]','fontsize',fsize);
ylabel('Northing [km]','fontsize',fsize);
title('Misfit Function','fontsize',fsize);
axis equal;     % WARNING: only for comparing x vs y
axis tight;

% plot solution - no mask
mogi2insar(xs(indx),ys(indy),zs,V,1,0);  
% MIGHT only show if imagesc is used in plot_model.m
subplot(2,1,1); hold on; plot(xs(indx),ys(indy),'kp','markersize',msize,'markerfacecolor','w');
subplot(2,1,2); hold on; plot(xs(indx),ys(indy),'kp','markersize',msize,'markerfacecolor','w');

% plot solution - with mask
syn_rngchg = mogi2insar(xs(indx),ys(indy),zs,V,1,1);  
% MIGHT only show if imagesc is used in plot_model.m
subplot(2,1,1); hold on; plot(xs(indx),ys(indy),'kp','markersize',msize,'markerfacecolor','w');
subplot(2,1,2); hold on; plot(xs(indx),ys(indy),'kp','markersize',msize,'markerfacecolor','w');

% plot residual field
residuals = obs_rngchg - syn_rngchg;
plot_model(residuals);

%==========================================================================
% fix xs and ys, search over zs and V

% implement V-zs search here


%==========================================================================
