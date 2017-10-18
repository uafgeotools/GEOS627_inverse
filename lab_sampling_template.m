%
% lab_sampling.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Example script for generating samples of a function using the rejection
% method (Tarantola 2005, Section 2.3.2).
%
% See notes for lab exercise lab_sampling.pdf
%
% This script is relevant to the epicenter homework problem that is based
% on Problem 7-1 in Tarantola (2005).
%
% calls plot_histo.m
%

clear, close all, clc

% limits for x to consider
% OUTSIDE THESE LIMITS WE ASSUME THAT p(x) = 0
xmin = -15; xmax = 12;

% define in-line function p(x)
% note that parameters like A1, A2, x1, x2 must be assigned before p()
x1 = -2; A1 = 2; sig1 = 2;
x2 =  4; A2 = 1; sig2 = 0.5;
ifun = 1;
switch ifun
    case 1
        p = @(x) ( A1*exp(-(x-x1).^2/(2*sig1^2)) );
    case 2
        p = @(x) ( A1*exp(-(x-x1).^2/(2*sig1^2)) + A2*exp(-(x-x2).^2/(2*sig2^2)) );
end

% KEY TECHNICAL POINT: f is a function, not a numerical array
% (note that x is not a stored array)
whos

% analytical curve
ncurve = 1000;
xcurve = linspace(xmin,xmax,ncurve);
pcurve = p(xcurve);

% generate samples
% KEY: what does rand do?
NTRY = 1e5;
xtry = xmin + (xmax-xmin)*rand(NTRY,1);

% sample the function
A = max([A1 A2]);           % note: only true for our choice of p(x)
%A = max(pcurve);            % (will work for any densely discretized p(x))
ptry = p(xtry) / A;         % SET A: values between 0 and 1
chance = rand(NTRY,1);      % SET B: values between 0 and 1

% plot
figure; nr=3; nc=2;
dxbin = 0.2;
dpbin = 0.05;    % 0.05, 0.1, 0.2
edgesx = [xmin:dxbin:xmax]; nex = length(edgesx);
edgesp = [0:dpbin:1];       nep = length(edgesp);

subplot(nr,nc,1); plot(xcurve,pcurve/A);
xlabel('x'); ylabel('p(x) / A'); title('(a)  p(x) / A'); axis([xmin xmax 0 1.2]);

subplot(nr,nc,2); plot_histo(xtry,edgesx); xlim([xmin xmax]);
xlabel('xtry'); title('(b)  xtry'); 

subplot(nr,nc,3); plot_histo(ptry,edgesp);
xlabel('p(xtry) / A'); title('(c)  p(xtry) / A'); 

subplot(nr,nc,4); plot_histo(chance,edgesp);
xlabel('ztry'); title('(d)  ztry'); 

% KEY COMMAND: compare pairs of test samples in sets A and B,
%              then accept or reject the test sample
ikeep = find(ptry > chance);
xkeep = xtry(ikeep);

subplot(nr,nc,5); plot_histo(xkeep,edgesx); xlim([xmin xmax]);
xlabel('xkeep'); title('(e)  xkeep');

% if p is a probability density and F is the misfit function, then
%    p(x) = exp(-F(x))
%    F(x) = -ln(p(x))
subplot(nr,nc,6); plot(xcurve,-log(pcurve));
axis([xmin xmax -1 1.1*max(-log(pcurve))]);
xlabel('x'); ylabel('F(x) = -ln(p(x))'); title('(f)  F(x)');

break  % first break

%==========================================================================
% ADDITIONAL PLOTTING TO UNDERSTAND THE PATTERN IN (c)

subplot(nr,nc,1); hold on;

ileftbin = nep - 1;
Pcut1 = edgesp(ileftbin);
Pcut2 = edgesp(ileftbin+1);
isub = find(and(pcurve/A >= Pcut1, pcurve/A < Pcut2 ));
plot(xcurve(isub),pcurve(isub)/A,'r.');
plot([xmin xmax],Pcut1*[1 1],'r--',[xmin xmax],Pcut2*[1 1],'r--');

ileftbin = round(nep*0.50);
Pcut1 = edgesp(ileftbin);
Pcut2 = edgesp(ileftbin+1);
isub = find(and(pcurve/A >= Pcut1, pcurve/A < Pcut2 ));
plot(xcurve(isub),pcurve(isub)/A,'k.');
plot([xmin xmax],Pcut1*[1 1],'k--',[xmin xmax],Pcut2*[1 1],'k--');

ileftbin = 1;
Pcut1 = edgesp(ileftbin);
Pcut2 = edgesp(ileftbin+1);
isub = find(and(pcurve/A >= Pcut1, pcurve/A < Pcut2 ));
plot(xcurve(isub),pcurve(isub)/A,'c.');
plot([xmin xmax],Pcut1*[1 1],'c--',[xmin xmax],Pcut2*[1 1],'c--');

%==========================================================================

break  % second break

% summary plot
figure; nr=2; nc=2; 
subplot(nr,nc,1); plot_histo(xtry,edgesx); xlim([xmin xmax]); xlabel('xtry');
subplot(nr,nc,2); plot_histo(ptry,edgesp); xlabel('p(xtry) / A');
subplot(nr,nc,3); plot_histo(xkeep,edgesx); xlim([xmin xmax]); xlabel('xkeep');
subplot(nr,nc,4); plot_histo(ptry(ikeep),edgesp); xlabel('p(xkeep) / A');

if bprint, fontsize(12);
    print(gcf,'-depsc',sprintf('%slab_sampling_hists3_ifun%i',pdir,ifun));
end

%==========================================================================
% PRACTICE WITH IN-LINE FUNCTIONS AND PLOTTING

break  % third break

% IMPLEMENT YOUR IN-LINE FUNCTION p HERE


% practice plotting this function
nx = 100;
ny = nx+1;      % avoid square arrays (for debugging)
n = nx*ny;
w = 2;
xvec = linspace(-w,w,nx);
yvec = linspace(-w,w,ny);
[X,Y] = meshgrid(xvec,yvec);
% OPTION A -- pass vector to p(), then plot as matrix
x = X(:);
y = Y(:);
pplot = p(x,y);
Pplot = reshape(pplot,ny,nx);
figure; pcolor(X,Y,Pplot); shading flat;
% OPTION B -- pass matrix to p(), then plot as matrix
Fplot = p(X,Y);
figure; pcolor(X,Y,Pplot); shading flat;
% OPTION C -- plot as vector
figure; scatter(x,y,4^2,pplot,'fill');

%==========================================================================
