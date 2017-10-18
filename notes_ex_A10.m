%
% notes_ex_A10.m
% Example A.10 of Aster et al., 2012 (2nd ed.)
% 

clc
clear
close all
format compact

disp('example matrix (A.49):');
A = [3 1 9 4 ; 2 1 7 3 ; 5 2 16 7]

disp('rref(A) (A.50):');
rref(A)

disp('basis vectors for null space (A.51):');
x1 = [-2 -3 1 0]'
x2 = [-1 -1 0 1]'

disp('a choice of b for A*x = b (A.52)');
b = [22 17 39]'         % given
%b = randi(99,3,1)       % most b are in the null space

rrefAB = rref([A b])

% choice of b leads to adding x0 to the solution: p = xn + x0
% (set x3 = x4 = 0)
x0 = [rrefAB(1:2,5) ; 0 ; 0]

% choose the "one particular solution" as in Aster
disp('choice of x3 and x4 lead to solution vector p (A.53):');
x3=1; x4=2;
% Eq A.53
xn = x3*x1 + x4*x2
disp('p = xn + x0:');
p = xn + x0

disp('add another null space vector to p (A.54):');
x3=2; x4=3;
%x3=0; x4=0;
xn = x3*x1 + x4*x2
disp('xn + p:');
xn + p

disp('check A*p - b:')
A*p - b
disp('check A*(xn+p) - b:')
A*(xn+p) - b

%==========================================================================
