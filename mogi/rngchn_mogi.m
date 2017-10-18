function [del_rng] = rngchn_mogi(xs,ys,zs,V,eing,ning,plook)
%RNGCHN_MOGI input mogi source and look vector, output look displacement
%
%  USEAGE: [del_rng] = rngchn_mogi(ys,xs,zs,del_v,ning,eing,plook);
%
%  INPUT: xs = local east coord of center of Mogi source (km)
%         ys = local north coord of center of Mogi source (km)
%         zs = elevation of Mogi source (km) (positive up)
%         V  = volume change of Mogi source (km^3)
%         eing = east coord's of points to calculate range change (x)
%         ning = north coord's of points to calculate range change (y)
%
%  OUTPUT: del_rng = range change at coordinates given in eing and ning.
%                    If eing and ning are vectors with the same dimensions,
%                    del_rng is a vector. If ning is a row vector and eing
%                    is a column vecor, del_rng is a matrix of deformation
%                    values...compliant with Feigle and Dupre's useage.
%
% CARL NOTE: The coding in this function could be shortened and improved.
%            (I doubt that two separate blocks are needed.)
%

[m,n] = size(ning);
[mm,nn] = size(eing);

% CARL: z is not allowed as an input coordinate, but I added it to the
% equations below to remind us that the displacement field DOES depend on
% elevation z within the elastic medium
z = 0;

%----coef for bulk modulus pressure <--> volume relation is below
%dsp_coef = (1000000*del_v*15)/(pi*16);

%----coef for McTigue's pressure <--> volume relation is below
% CARL: presumably this is Eq 1: (1/pi)*(1-nu)*del_v,
%       where nu = 0.25 is for a Poisson medium and
%       1e6 converts from km to mm
dsp_coef = (1e6*V*3)/(pi*4);

if(mm==1 && n==1)
    %disp('Calculating a matrix of rngchg values')
    del_rng = zeros(m,nn);
    %del_d = del_rng;
    %del_f = del_rng;
    tmp_n = del_rng;
    tmp_e = del_rng;
    for i_loop = 1:m
        tmp_e(i_loop,:) = eing;
    end
    for i_loop = 1:nn
        tmp_n(:,i_loop) = ning;
    end
    % horizontal distance
    d_mat = sqrt((tmp_n-ys).^2 + (tmp_e-xs).^2);
    % denominator of displacement field for mogi source
    tmp_hyp = ((d_mat.^2 + (z-zs).^2).^1.5);
    % Uh, horizontal displacement
    uh = dsp_coef*d_mat./tmp_hyp;
    % Uz
    uz = -dsp_coef*zs./tmp_hyp;
    % azimuthal angle
    azim = atan2((tmp_e-xs),(tmp_n-ys));
    % Ux and Uy
    e_disp = sin(azim).*uh;
    n_disp = cos(azim).*uh;
    % project displacement field onto look vector
    for i_loop = 1:nn
        del_rng(:,i_loop) = ...
            [e_disp(:,i_loop) n_disp(:,i_loop) uz(:,i_loop)]*plook';
    end
    % flip sign projected distance (or flip plook)
    del_rng = -1.0*del_rng;
    
elseif ((mm==1 && m==1) || (n==1 && nn==1))
    if (n~=nn)
        error('Coord vectors not equal length!')
    end
    %disp('Calculating a vector of rngchng values')
    % horizontal distance
    d_mat = sqrt((ning-ys).^2 + (eing-xs).^2);
    % denominator of displacement field for mogi source
    tmp_hyp = ((d_mat.^2 + (z-zs).^2).^1.5);
    % Uh, horizontal displacement
    uh = dsp_coef*d_mat./tmp_hyp;
    % Uz
    uz = -dsp_coef*zs./tmp_hyp;
    % azimuthal angle
    azim = atan2((eing-xs),(ning-ys));
    % Ux and Uy
    e_disp = sin(azim).*uh;
    n_disp = cos(azim).*uh;
    % project displacement field onto look vector
    del_rng = [e_disp n_disp uz]*plook';
    % flip sign of projected distance (or flip plook)
    del_rng = -1.0*del_rng;
    
else
    error('Coord vectors make no sense!')
end
