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

#   Source:
#   W. Li and J. Swetits,
#   "A Newton method for convex regression, data smoothing and
#   quadratic programming with bounded constraints",
#   SIAM J. Optimization 3 (3) pp 466-488, 1993.

#   SIF input: Nick Gould, August 1994.

#   classification QLR2-AN-V-V

Translated to scilab from AMPL by Sharvani Laxmi Somayaji and Pranshu Malhotra as a part of FOSSEE internship, 2021
*/

/*
-------test case outputs-------
    
For N=10, in scilab
fval  = 0.0219117
 output  = 

  Iterations = 27
  Cpu_Time = 0.066
  Objective_Evaluation = 37
  Dual_Infeasibility = 1.110D-16
  Message = "Optimal Solution Found"

For same N, on NEOS server:
objective = 2.1911074660776464e-02 (scaled and unscaled)

For N=1000, in scilab
fval  = 4.9316211
 output  = 

  Iterations = 31
  Cpu_Time = 130.769
  Objective_Evaluation = 32
  Dual_Infeasibility = 1.043D-12
  Message = "Optimal Solution Found"
For same N, on NEOS server:
objective = 4.9311300941879574e+00 (scaled and unscaled)
*/

funcprot(0);
//N = 1000;
N=10; 
    
K = 2;

B = zeros(1,K+1);
B(1) = 1;
for i = 2:K+1
    B(i) = (i-1) * B(i-1);
end

C = zeros(1,K+1);
C(1) = 1;
i = 2:(K+1);
C(i) = (-1).^(i-1).*B(K+1)./(B(i).*B(K+1-i+1));

x0 = zeros(1,N+K);

//linear inequality constraints
j =1:N;
i = 0:K;
A = zeros(N,(N+K)); 
b=zeros(N,1);
for j=1:N
    A(j,j+K-i) = -C(i+1);
end

function y=f(x)
    //N = 1000;
    N=10;
    
    K = 2;
    T = zeros(1,N+K);
    i = 1:(N+K);
    T(i) = (i-1)/(N+K-1);
    pi = 3.1415;
    term1 = sum(-(cos(pi.*T(i)) + 0.1*sin(i)).*x(i));
    term2 = sum(0.5.*( cos(pi.*T(i)) + 0.1.*sin(i) ).^2);
    term3 = sum(0.5.*x(i).^2);

    y = term1 + term2 + term3;
endfunction

function g = fGrad(x)
    //N = 1000;
    N=10;
    K = 2;
    T = zeros(1,N+K);
    i = 1:(N+K);
    T(i) = (i-1)/(N+K-1);
    pi = 3.1415;
    
    g = zeros(N+K, 1);
    g(i,1) = -(cos(pi.*T(i)) + 0.1*sin(i)) + x(i);

endfunction

options = struct("MaxIter", [1000000], "CpuTime", [60000], "GradObj", fGrad, "Hessian","off","GradCon","off","HessianApproximation", [1] );
[x,fval,exitflag,output] =fot_fmincon(f, x0,A,b,[],[],[],[],[],options)
