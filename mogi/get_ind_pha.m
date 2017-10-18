function  [eing_vec,ning_vec,vec_ind] = get_ind_pha(pha,ning,eing)
%convert a 2-D image into three 1-D arrays
% USEAGE: [eing_vec,ning_vec,vec_ind]=get_ind_pha(pha,ning,eing);
%
%  INPUT:  pha = phase matrix
%          eing = easting vector
%          ning = (reversed) northing vector
%
%  OUTPUT: eing_vec = vector of eastings corres. to pha in vector form
%          ning_vec = vector of northings corresp. to pha  in vector form
%          vec_ind  = index of vectors to be run through rngchn.f

[m,n] = size(pha);
mm = length(ning);
nn = length(eing);
if(mm ~= m || nn ~= n)
  error('I have to leave GET_IND')
end

pha_vec = reshape(pha,m*n,1);
vec_ind = find(pha_vec ~= -9999.0);

[nck,mck] = size(ning);
if (nck == 1) 
  ning = ning';
end
tmp_mat1 = zeros(m,n);
for i=1:n
  tmp_mat(:,i)=ning; 
end
tmp_vec = reshape(tmp_mat,m*n,1);
ning_vec = tmp_vec(vec_ind);

[nck,mck] = size(eing);
if (mck == 1) 
  eing=eing';
end
for i=1:m
  tmp_mat(i,:) = eing;
end
tmp_vec = reshape(tmp_mat,m*n,1);
eing_vec = tmp_vec(vec_ind);
