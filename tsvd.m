function [f_r, rss, f_r_ss] = tsvd(g, X, rvec)
%TSVD regularization using truncated singular value decomposition
%
% INPUT
%   g       n x 1 data vector
%   X       n x p design matrix
%   rvec    r x 1 vector of truncation parameters (integers between 1 and p)
%
% OUTPUT
%   f_r     p x r matrix of TSVD model vectors
%   rss     r x 1 vector of residual sum of squares
%   f_r_ss  r x 1 vector of norm-squared of each f_r vector
%
% Given a vector g, a design matrix X, and a truncation parameter r, 
%       [f_r, rss, f_r_ss] = tsvd(g, X, r) 
% this returns the truncated SVD estimate of the vector f in the linear
% regression model
%        g = X*f + noise
% 
% If r is a vector of truncation parameters, then the ith column
% f_r(:,i) is the truncated SVD estimate for the truncation
% parameter r(i); the ith elements of rss and f_r_ss are the
% associated residual sum of squares and estimate sum of squares.
% 
% Adapted from TSVD routine in Per Christian Hansen's Regularization Toolbox. 
%

% size of inputs (n is number of data; p is number of parameters)
[n, p]      = size(X);
q           = min(n, p);
nr          = length(rvec);

% Possible choice of truncation parameter?
if ( min(rvec) < 0 || max(rvec) > q )
    rvec, q
    error('Impossible truncation parameter r.')
end

% initialize outputs
f_r         = zeros(p, nr);             % set of r models
rss         = zeros(nr, 1);             % RSS for each model
f_r_ss      = zeros(nr, 1);             % norm of each model

% compute SVD of X
[U, S, V]   = svd(X, 0);  
s           = diag(S);                  % vector of singular values

% 'Fourier' coefficients (fc) in expansion of solution in terms of right singular vectors
% note: these are also referred to as Picard ratios
beta        = U(:, 1:q)'*g;             % note data g
fc          = beta ./ s;

% treat each truncation parameter separately
for j = 1:nr
	k         = rvec(j);                % current truncation parameter
	f_r(:, j) = V(:, 1:k) * fc(1:k);    % truncated SVD estimated model vector
	f_r_ss(j) = sum(fc(1:k).^2);        % the squared norm of f_r
	rss(j)    = sum(beta(k+1:q).^2);    % residual sum of squares
end

% in overdetermined case, add rss of least-squares problem
if (n > p)
    rss = rss + sum((g - U(:, 1:q)*beta).^2);   % note data g
end

%==========================================================================
  