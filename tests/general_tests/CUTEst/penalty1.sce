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

#   Source:  Problem 23 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley #181 (p. 79).

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

Translated to scilab from AMPL by Sharvani Laxmi Somayaji and Pranshu Malhotra as a part of FOSSEE internship, 2021
*/

/*
---------test case outputs--------
**Note that the outputs will depend on various factors such as solver, the problem type, tolerance etc**

N=10; In scilab:
fval = 0.0001305;
Cpu_Time: 0.076
 output  = 

  Iterations = 14
  Cpu_Time = 0.076
  Objective_Evaluation = 22
  Dual_Infeasibility = 0.0000189
  Message = "Optimal Solution Found"
  
  for same n in NEOS server, fval = 4.6053629348787437e-07 (scaled) 
  for same n in NEOS server, fval = 7.0876536396749192e-05 (unscaled) 
  f = 7.08765e-05

N=100; In scilab:
fval = 0.0635000;
Cpu_Time: 0.008
 output  = 

  Iterations = 1
  Cpu_Time = 0.008
  Objective_Evaluation = 2
  Dual_Infeasibility = 0.00002
  Message = "Optimal Solution Found"

  for same n in NEOS server, fval = 7.4530367665711525e-10 (scaled) 
  for same n in NEOS server, fval = 1.0086932506988202e-03 (unscaled) 
  f = 0.00100869

N=1000; In scilab:
  fval  = 2.1897404
 output  = 

  Iterations = 32
  Cpu_Time = 0.188
  Objective_Evaluation = 40
  Dual_Infeasibility = 0.4246078
  Message = "Optimal Solution Found"

For same n in gamsworld site list, fval = 0.00968618

for same n in NEOS server, fval or objective = 6.4394978143979231e-08 (scaled)
for same n in NEOS server, fval = 6.4394978143979227e+00(unscaled)
f = 6.4395
*/

funcprot(0);
//N=1000;
N=100;
M = N+1;

x0 = zeros(1,N);
i = 1 : N;
x0(1,i) = i;
a = 10^-5;

function y=f(x)
    //N=1000;
    N=100;
    M = N+1;
    a = 10^-5;
    j = 1:N;
    term1 =( sum((x(j)).^2)  - 0.25 ).^2; 
    term2 =sum(a .*( (x(j) - 1).^2) );
  
    y = term1+ term2;
endfunction

function g=fGrad(x)
    //N=1000;
    N=100;
    M = N+1;
    a = 10^-5;
    j = 1:N;
    sqsum = sum( (x(j)).^2 );
    
    g = zeros(1,N);
    g(1,j) = (4.*x(j)).*(sqsum - 0.25) + 2.*a.*(x(j)-1) ;  
endfunction

options = struct("MaxIter", [1000000], "CpuTime", [600000], "GradObj", fGrad, "Hessian","off","GradCon","off","HessianApproximation", [1] );
[x,fval,exitflag,output] =fot_fmincon(f, x0,[],[],[],[],[],[],[],options)
