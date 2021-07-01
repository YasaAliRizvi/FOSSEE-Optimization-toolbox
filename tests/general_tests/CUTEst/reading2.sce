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
#   S. Lyle and N.K. Nichols,
#   "Numerical Methods for Optimal Control Problems with State Constraints",
#   Numerical Analysis Report 8/91, Dept of Mathematics, 
#   University of Reading, UK.

#   SIF input: Nick Gould, July 1991.

#   classification LLR2-MN-V-V

Translated to scilab from AMPL by Yasa Ali Rizvi as a part of FOSSEE internship, 2021
*/

/*
-------test case outputs-------
//  Bestknown Objective = -0.00996705 (= -0.0125523; N=300)

//N = 100; //fval =  -0.0124710
N = 300; //fval =   -0.0124547

For N=50 in scilab, 
 fval  = -0.0123436
 output  = 

  Iterations = 357
  Cpu_Time = 1.428
  Objective_Evaluation = 360
  Dual_Infeasibility = 4.996D-16
  Message = "Optimal Solution Found"
for same N value on NEOS server, 
    objective = -1.2353115577435170e-02
*/
funcprot(0);
//N = 100; 
N = 50; 
A = 0.07716;
H = 1/N;
pi = 3.1415;

x0 = zeros(3*N+3, 1);

//lower and upper bounds
lb = zeros(1, 3*N+3);;
ub = zeros(1, 3*N+3);;

i=1:(N+1);
lb(i) = -%inf;
ub(i) = %inf;

i=(N+2):(2*N+2);
lb(i) = -0.125;
ub(i) = 0.125;

i=(2*N+3):(3*N+3);
lb(i) = -1.0;
ub(i) = 1.0;

//Linear equality constraints
Aeq = [];
beq = zeros(2*N+2,1);

//constraint 1
for i=1:N
    Aeq(i,i) = -1/H;
    Aeq(i,i+1) = 1/H;
    Aeq(i,N+1+i) = -0.5;
    Aeq(i,N+1+i+1) = -0.5;
end

//constraint 2
for i=N+1:2*N
    Aeq(i,i+1) = -1/H;
    Aeq(i,i+2) = 1/H;
    Aeq(i,N+1+i+1) = -0.5;
    Aeq(i,N+1+i+2) = -0.5;
end

//constraint 3
Aeq(2*N+1, 1) = 1;

//constraint 4
Aeq(2*N+2, N+2) = 1;

function y=f(x)
    N=50;
    H = 1/N;
    pi = 3.1415;
    sum = 0;
    for i=1:N
        sum = sum + (-0.5*H*cos(2*pi*i*H)*x(i+1) - ..
        0.5*H*cos(2*pi*(i-1)*H)*x(i+1) + H*(x(2*N+2+i+1)+x(2*N+2+i))/(8*pi^2));
    end
    y = sum;
endfunction

function y=fGrad(x)
    N=50;
    A = 0.07716;
    H = 1/N;
    pi = 3.1415;
    y(1,1) = 0;
    for i=1:N
        y(1,i+1) = -0.5*H*cos(2*pi*i*H) - 0.5*H*cos(2*pi*(i-1)*H);
    end
    y(1,2*N+3) = H/(8*pi^2);
    y(1,3*N+3) = H/(8*pi^2);
    for i=2*N+4:3*N+2
        y(1,i) = H/(4*pi^2);
    end
endfunction

options = struct("MaxIter", [3000], "CpuTime", [600], "HessianApproximation", [1], "GradObj", fGrad, "Hessian","off","GradCon","off" );
[x,fval] = fot_fmincon(f, x0, [], [], Aeq, beq, lb, ub, options)
