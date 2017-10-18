function [real_matrix]=read_real_be(fname,num_cols,num_rows);
%
%  USEAGE: [real_matrix]=read_real_be(fname,num_rows,num_cols)
%
%  INPUT: fname = name of coherence file
%         num_rows = number of rows in coherence matrix
%         num_cols = number of columns " "
%
%  OUTPUT: real_matrix = matrix of floating point numbers
%

%-----First make sure the file exists:
if ~exist(fname) 
  mess=['File: "' fname '" Does not exist!'];
  disp(mess)
  error('I"m getting out of here!')
end

disp(['Reading values from ',fname])
 
%-----Open the file and read everything:
fid =fopen(fname,'r', 'b');
[coh,count]=fread(fid,'float');
st=fclose(fid);
if st == -1
  mess=['Now why can"t I close this file?'];
  disp(mess)
end
 
%-----check below to see if the dimensions fed by num_cols and num_rows
%-----are consistent with the number of phase values read from the file:
icheck=num_rows*num_cols;
if icheck ~= count
  mess=['Dimensions given are different from those read'];
  disp(mess)
  error('I have to go!')
end

%-----find coherence values less than coh_cutoff and change to NaN

%-----crop the matrix as predetermined  and pass the values back
real_matrix=reshape(coh,num_rows,num_cols);
real_matrix=real_matrix';

