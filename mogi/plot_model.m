function [xvec,yvec] = plot_model(D,bonefig)
%PLOT_MODEL plot an interferogram both as an unwrapped image and wrapped phase image

SAMPLE = 1100;
LINE = 980;
POSTING = 40.0;
%HALF_WAVE = 28.3;

xvec = [1:POSTING:SAMPLE*POSTING]/1000;
yvec = [1:POSTING:LINE*POSTING]/1000;
yvec = fliplr(yvec);

F1 = D;
F2 = wrap(D/10)/5.66*4*pi;

bimagesc = false;       % plot with imagesc (=true) or pcolor (=false)
nanclr = 0.6*[1 1 1];   % NaN color (imagesc=0 only)

% =true for one figure; =false for two figures
if nargin==1, bonefig = true; end
if bonefig==true; figure; end

for kk=1:2

    if kk==1
        F = D;
        %tlab = 'Line of Sight Motion [mm]';
        tlab = 'Displacement in look direction [mm]';
        clims = [-30 30];
    else
        F = wrap(D/10)/5.66*4*pi;
        tlab = 'Interferogram phase [rad]';
        clims = [-1 1]*max(abs(F(:)));
    end
    
    if bonefig, subplot(2,1,kk); else figure; end
    if bimagesc
        imagesc(xvec,yvec,F,clims);
        set(gca,'ydir','normal');
    else
        pcolor(xvec,yvec,F); shading flat;
        %set(gca,'ydir','reverse');
        caxis(clims);

        % silly Matlab commands to make sure that the colors that are shown in
        % the figure are actually printed to files
        set(gcf,'Color','white');
        set(gca,'Color',nanclr);
        set(gcf,'InvertHardCopy','off');
    end
    colorbar   
    axis equal; axis tight
    set(gca,'FontSize',12);
    xlabel('Easting [km]','FontSize',14);
    ylabel('Northing [km]','FontSize',14);
    title(tlab,'FontSize',14);

end
