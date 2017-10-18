function Cd = covC(id,parms)
%COVC evaluate covariance function C(d) at an array of distances d
%
% INPUT
%   id      array of distance values (scalar, vector, matrix, etc)
%   parms:
%   icov    type of covariance function (=1 Gaussian; =2 exponential)
%   iL      length scale (same units as id)
%   sigma   amplitude factor
%   nu      OPTIONAL: parameter for Matern covariance (icov=4 only)
%
% OPTIONS FOR SPECIFYING LENGTH SCALE
%   (1) iL and id are indices for a spatial grid
%   (2) iL and id are actual lengths for a spatial grid
%
% EXAMPLE: d = linspace(0,100); Cd = covC(d,{1,30,3}); figure; plot(d,Cd); grid on;
% 
% For context for the first theree covariance functions,
% see Tarantola (2005), Section 5.3.3 (p. 113).
% To match his convention for L, set LFACTOR = 1.
% My preference is to use LFACTOR = 2, since this provides invariance
% of C(L) for covariance functions within the Matern family (vary nu).
% For the circular covariance function, our L is Tarantola's D.
%
% Carl Tape, 15-Jan-2010
%

% get input parameters
nparm = length(parms);
icov  = parms{1};
iL    = parms{2};
sigma = parms{3};
nu    = [];
if nparm==4, nu = parms{4}; end

if ~any(icov==[1 2 3 4]), error('icov must be 1,2,3,4'); end
if and(isempty(nu),icov==4)
   error('must provide nu as input for Matern covariance (icov=4)');
end

% see text above
LFACTOR = 2;    % our preference
%LFACTOR = 1;    % formulas in Tarantola 2005
%disp(sprintf('covC.m: LFACTOR = %.1f',LFACTOR));

% compute covariance for input distances (id)
switch icov
    case 1
        % Gaussian covariance
        % --> The factor of 2 in (2*iL^2) leads to smoother models
        iL = iL / LFACTOR;
        Cd = sigma^2 * exp(-id.^2 / (2*iL^2) );
        
    case 2
        % exponential covariance
        iL = iL / LFACTOR;
        Cd = sigma^2 * exp(-id / iL );
        
    case 3
        % circular covariance
        % here iL represents the diameter of the two intersecting discs
        icirc = find(id <= iL);
        Cd = zeros(size(id));
        beta = 2*asin(id(icirc) / iL);
        Cd(icirc) = sigma^2 * (1 - (beta + sin(beta))/pi );
        
    case 4
        % Matern covariance
        % http://en.wikipedia.org/wiki/Mat%C3%A9rn_covariance_function
        % note this uses the built-in functions gamma and besselk
        iL = iL / LFACTOR;
        b = besselk(nu, sqrt(2*nu)*id/iL);
        Cd = sigma^2 * (1/(gamma(nu)*2^(nu-1))) * (sqrt(2*nu)*id/iL).^nu .* b;
end
