% Generates n copies of a complex 1-D Gaussian random field 
% from its spectral distribution matrix

% Input: 
% k = vector of k values
% C = spectral distribution vector (must be a row vector). 
% n = number of copies to generate

% Output
% phi = complex Gaussian random field with covariance function 
%       mhifft(C)


function phi = grf1(k, C, n)

Nx = max(size(k));
dk = k(2) - k(1);
Periodx = 2*pi/dk;

Cmtx = repmat(C,[n,1]);  % Generate copies of the spectral matrix
A = randn([n,Nx]);   % N(0,1) random variables
B = randn([n,Nx]);

% Generate random field 
phi = sqrt(Periodx*Cmtx/2).*(A + 1i*B);


