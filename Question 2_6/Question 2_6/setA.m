function A = setA(vecA)
% input: vector of unique elements of matrix A
% output: full matrix A

%param = [alpha_pipi; rho; trinche_y; trinche_pi; alpha_my; alpha_mpi; alpha_mr]

%param = [alpha_pipi; alpha_ry; alpha_rpi; alpha_my; alpha_mpi; alpha_mr]

A = [1.0        0.0        0.0        0.0; ...
     -vecA(1)   1.0        0.0        0.0; ...
     -vecA(2)   -vecA(3)   1.0        0.0; ...
     -vecA(4)   -vecA(5)   -vecA(6)   1.0];
end


