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

#   Source: problem P4 in
#   W.W. Hager,
#   "Multiplier Methods for Nonlinear Optimal Control",
#   SIAM J. on Numercal Analysis 27(4): 1061-1080, 1990.

#   SIF input: Ph. Toint, April 1991.

#   classification OLR2-AN-V-V

Translated to scilab from AMPL by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
*/

/*
--------test case outputs-----------
For N=50, in scilab

 fval  = 2.7998374
output  = 

  Iterations = 13
  Cpu_Time = 0.061
  Objective_Evaluation = 21
  Dual_Infeasibility = 1.580D-09
  Message = "Optimal Solution Found"

For same N, on NEOS server:
objective = 2.7998109577397630e+00

For N=500, in scilab
fval  = 2.7945138
output  = 

  Iterations = 10
  Cpu_Time = 8.621
  Objective_Evaluation = 11
  Dual_Infeasibility = 0.0000004
  Message = "Optimal Solution Found"

For same N, on NEOS server:
fval = 2.7945134554215989e+00

*/

//n = 5000;

//n=500;
n=50;

h = 1/n;
i=0:n;
t= zeros(n+1);z= zeros(n+1);
t(i+1) = i*h; z(i+1) = exp(-2*t(i+1));
i = 1:2;
a(i) = -0.5*z(i);
b(i) = a(i).*(t(i)+0.5);
c(i)= a(i).*(t(i).^2+t(i)+0.5);

scda = (a(2)-a(1))/2;
scdb = (b(2)-b(1))*n;
scdc = (c(2)-c(1))*n*n*0.5;

e = exp(1);
xx0 = (1+3*e)/(2-2*e);

x0 = zeros(1,(n+1+n));

xx(1) = xx0;
x0(1) = xx(1);

//Linear equality constraints
Aeq = zeros( n+1, 2*n+1 ); beq = zeros(n+1,1);
Aeq(1,1) = 1; beq(1,1) = xx0;
j = 2:(n+1);
for j = 2:(n+1)
    Aeq(j,j-1) = -n;
    Aeq(j,j) = n-1;
    Aeq(j,n+j) = -exp(t(j));
end

//lower and upper bounds
ub = zeros(1,2*n+1);
ub(1,1:(2*n+1))=%inf;
ub(1,(n+2):(2*n+1)) = 1;
    
function y=f(x)
    //n = 5000;
    //n=500;
    n=50;
    
    h = 1/n;
    i=0:n;
    t= zeros(n+1);z= zeros(n+1);
    t(i+1) = i*h; z(i+1) = exp(-2.*t(i+1));
    i = 1:2;
    a(i) = -0.5*z(i);
    b(i) = a(i).*(t(i)+0.5);
    c(i)= a(i).*(t(i).^2+t(i)+0.5);
    
    scda = (a(2)-a(1))/2;
    scdb = (b(2)-b(1))*n;
    scdc = (c(2)-c(1))*n*n*0.5;
    
    xx = x(1: (n+1));
    xu = x((n+1+1):(n+1+n));
    
    term2 = sum((xu.^2)*h*0.5);
    i=2:(n+1)
    term1 = sum( (scda.*z(i-1).*xx(i).^2) + scdb.*z(i-1).*xx(i).*(xx(i-1)-xx(i)) + scdc.*z(i-1).*(xx(i-1)-xx(i)).^2 );
    y = term1 + term2;
endfunction

function g = fGrad(x)
    //n = 5000;
    //n=500;
    n=50;
    
    //constants which don't need changing
    h = 1/n;
    i=0:n;
    t= zeros(n+1);z= zeros(n+1);
    t(i+1) = i*h; z(i+1) = exp(-2.*t(i+1));
    i = 1:2;
    a(i) = -0.5*z(i);
    b(i) = a(i).*(t(i)+0.5);
    c(i)= a(i).*(t(i).^2+t(i)+0.5);
    scda = (a(2)-a(1))/2;
    scdb = (b(2)-b(1))*n;
    scdc = (c(2)-c(1))*n*n*0.5;
    
    xx = x(1: (n+1));
    xu = x((n+1+1):(n+1+n));
    usum = sum((xu.^2)*h*0.5);
    y = 0;
    i=2:(n+1)
        
    gu = zeros(n); i = 1:n;
    gu(i) = 2.*(xu(i)).*(0.5*h);
    gx = zeros(n+1); i = 2:n;
    gx(i) = 2.*xx(i).*((scda.*z(i-1))-(scdb.*z(i-1))+(scdc.*z(i-1))) + 2.*xx(i).*((scdc.*z(i))) + xx(i-1).*((scdb.*z(i-1))-2.*(scdc.*z(i-1))) + xx(i+1).*((scdb.*z(i-1))-2.*(scdc.*z(i))) ;
    //for i=1
    i=1;
    gx(i) = 2.*xx(i).*((scdc.*z(i))) + xx(i+1).*((scdb.*z(i))-2.*(scdc.*z(i))) ;
    //for i=n+1
    i = (n+1);
    gx(n+1) = 2.*xx(i).*((scda.*z(i-1))-(scdb.*z(i-1))+(scdc.*z(i-1))) + xx(i-1).*((scdb.*z(i-1))-2.*(scdc.*z(i-1))) ;
    
    g = zeros(2*n+1,1);
    i=1:(n+1);
    g(i,1) = gx(i);
    i = 1:n;
    g((n+1+i),1) = gu(i);
endfunction

options = struct("MaxIter", [1000000], "CpuTime", [60000], "GradObj", fGrad, "Hessian","off","GradCon","off","HessianApproximation", [1] );
[x,fval,exitflag,output] =fot_fmincon(f, x0,[],[],Aeq,beq,[],ub,[],options)
