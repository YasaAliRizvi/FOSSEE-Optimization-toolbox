/*
# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source: problem 3 in
#   M. J. D. Powell,
#   "Log barrier methods for semi-infinite programming calculations"
#   Numerical Analysis Report DAMTP 1992/NA11, U. of Cambridge, UK.

#   SIF input: A. R. Conn and Nick Gould, August 1993

#   classification LLR2-AN-4-V

Translated to scilab from AMPL by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
*/

/*
-------------test case outputs--------------

For M=100, in scilab
CPU time: 0.168
fval = 0.4964553
For same M, on NEOS server:
fval = 4.9645445145073652e-01


For M=1000, in scilab
CPU time: 0.395
fval = 0.5334123
For same M, on NEOS server:
fval = 5.3341145259338096e-01

For M=10000, in scilab
CPU time: 2.988
fval = 0.5356508
For same M, on NEOS server:
fval = 5.3565076485080632e-01

*/

//M=1000;
M=10000;
STEP = 8/M;
xi = zeros(M/2);
for j=1:(M/2)
    if (1 <= j && j<= M/8)
        xi(j) = 0;
    elseif (M/8+1 <= j && j<= M/4)
        xi(j) = (j-1)*STEP-1;
    elseif (M/4+1 <= j && j<= 3*M/8)
        xi(j) = 1;
    elseif (3*M/8+1 <= j && j<= M/2)
        xi(j) = (j-1)*STEP-3;
    end
end
eta = zeros(M/2);
for j=1:(M/2)
    if (1 <= j && j<= M/8)
        eta(j) = (j-1)*STEP;
    elseif (M/8+1 <= j && j<= M/4)
        eta(j) = 1;
    elseif (M/4+1 <= j && j<= 3*M/8)
        eta(j) = (j-1)*STEP-2;
    elseif (3*M/8+1 <= j && j<= M/2)
        eta(j) = 0;
    end
end

x_init = zeros(4, 1);
x_init = [-0.1;0;0;1.2];
x0 = x_init;

//linear inequality constraints
//Constraint 1
j=1:(M/2);
A1 = zeros(M/2, 4); b1=zeros(M/2,1);
A1(j,1) = -1; 
A1(j,2) = -xi(j); 
A1(j,3) = -eta(j); 
A1(j,4) = -1;
b1(j) = -xi(j).^2.*eta(j);
//Constraint 2
A2 = zeros(M/2, 4); b2=zeros(M/2,1);
A2(j,1) = 1; 
A2(j,2) = xi(j); 
A2(j,3) = eta(j); 
A2(j,4) = 0;
b2(j) = xi(j).^2.*eta(j);

A = [A1;A2]; b = [b1;b2];

function y = f(x)
	y = x(4);
endfunction

function g = fGrad(x)
    g = zeros(4,1);
    g(4,1) = 1;
endfunction

options = struct("MaxIter", [1000000], "CpuTime", [600000], "GradObj", fGrad, "Hessian","off","GradCon","off","HessianApproximation", [1] );
[x,fval,exitflag,output] =fot_fmincon(f, x0,A,b,[],[],[],[],[],options)
