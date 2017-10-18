%
% lab_svdgeom.m
% Carl Tape, GEOS 627, Inverse Problems and Parameter Estimation
%
% Geometry of the singular value decomposition
%

clc, clear, close all
 
% pick a matrix
G = randi([-2 2],2,2)

% singular value decomposition
[U,S,V] = svd(G)
Gcheck = U*S*V'

% singular values
s1 = S(1,1);
s2 = S(2,2);
% basis vectors
v1 = V(:,1);
v2 = V(:,2);
u1 = U(:,1);
u2 = U(:,2);

% unit circle for plotting
R = 1;
n = 100;
t = linspace(0,2*pi,n);
x = R*cos(t);
y = R*sin(t);
Cx = [x ; y];

% unit circles representing basis U and basis V
Ux = Cx;
Vx = Cx;

% standard basis
e1 = [1 0]';
e2 = [0 1]';

ax0 = 3*[-1 1 -1 1];
figure; nr=2; nc=2; fsize = 12;

subplot(nr,nc,1); hold on; grid on;
plot(Vx(1,:),Vx(2,:),'b-'); axis equal, axis(ax0);
plot([0 v1(1)],[0 v1(2)],'b',[0 v2(1)],[0 v2(2)],'b','linewidth',2);

GVx = G*Vx;     % ellipse for plotting
su1 = s1*u1;
su2 = s2*u2;
Ge1 = G*e1;
Ge2 = G*e2;

subplot(nr,nc,2); hold on; grid on;
plot(GVx(1,:),GVx(2,:),'r'); axis equal, axis(ax0);
plot(Ux(1,:),Ux(2,:),'r--');
plot([0 su1(1)],[0 su1(2)],'r',[0 su2(1)],[0 su2(2)],'r','linewidth',2);

% ugly code for positioning the text labels
f = 0.5;
subplot(nr,nc,1); 
[th,r] = cart2pol(v1(1),v1(2)); [v1p(1),v1p(2)] = pol2cart(th,r+f);
[th,r] = cart2pol(v2(1),v2(2)); [v2p(1),v2p(2)] = pol2cart(th,r+f);
text(v1p(1),v1p(2),'v_1','fontsize',fsize,'horizontalalignment','center');
text(v2p(1),v2p(2),'v_2','fontsize',fsize,'horizontalalignment','center');
subplot(nr,nc,2); 
[th,r] = cart2pol(su1(1),su1(2)); [su1p(1),su1p(2)] = pol2cart(th,r+f);
[th,r] = cart2pol(su2(1),su2(2)); [su2p(1),su2p(2)] = pol2cart(th,r+f);
text(su1p(1),su1p(2),'s_1 u_1','fontsize',fsize,'horizontalalignment','center');
text(su2p(1),su2p(2),'s_2 u_2','fontsize',fsize,'horizontalalignment','center');

% START CODE HERE FOR U --> GDAGGER --> V


%==========================================================================
