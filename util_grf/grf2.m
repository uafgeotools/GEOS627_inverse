function [phi,A,B] = grf2(k,m,C,n,A,B)
%GRF2 generate n copies of a complex Gaussian random field from its spectral distribution matrix
%
% INPUT
%   k      vector of k values
%   m       vector of m values
%   C       spectral distribution matrix. Assumed to have been generated
%           from meshgrid(K,M), so that it is indexed as C(m,k)
%   n       number of copies to generate
%   A,B     OPTIONAL: gaussian random vectors
%
% OUTPUT
%   phi     complex Gaussian random field with covariance function mhifft(C)
%   A,B     gaussian random vectors
%
% Two ways to use this function. Consider that you have two different
% covariance functions, C1 and C2.
%   option 1. use different random vectors to generate the GRFs
%      phi1 = grf2(k,m,C1,n)
%      phi2 = grf2(k,m,C2,n)
%   option 2. use the same random vectors to generate the GRFs
%      [phi1,A,B] = grf2(k,m,C1,n)
%      phi2 = grf2(k,m,C2,n,A,B)
%

Nx = max(size(k));
Ny = max(size(m));
dk = k(2) - k(1);
dm = m(2) - m(1);
Periodx = 2*pi/dk;
Periody = 2*pi/dm;

Cmtx = repmat(C,[1,1,n]);  % Generate copies of the spectral matrix

if nargin==6
    disp('grf2mod.m: using the provided Gaussian random vectors (A and B)'); 
else
    A = randn([Ny,Nx,n]);   % N(0,1) random variables
    B = randn([Ny,Nx,n]);
    %if nargin==5, phi = []; return; end
end
    
% generate Gaussian random field (covariance prescribed by Cmtx)
phi = sqrt(Periodx*Periody*Cmtx/2).*(A + 1i*B);
