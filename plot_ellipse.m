%set the number of points on the ellipse to generate and plot
function [x,y] = plot_ellipse(DELTA2,C,m)
n=1000;
%construct a vector of n equally-spaced angles from (0,2*pi)
theta=linspace(0,2*pi,n)';
%corresponding unit vector
xhat=[cos(theta),sin(theta)];
Cinv=inv(C);
%preallocate output array
r=zeros(n,2);
for i=1:n
%store each (x,y) pair on the confidence ellipse
%in the corresponding row of r
r(i,:)=sqrt(DELTA2/(xhat(i,:)*Cinv*xhat(i,:)'))*xhat(i,:);
end
x = m(1)+r(:,1);
y = m(2)+r(:,2);
plot(x,y);

% CHT: what is the justification for this command?
% clearly it has a big impact on what appears to be correlated or not
%axis equal
