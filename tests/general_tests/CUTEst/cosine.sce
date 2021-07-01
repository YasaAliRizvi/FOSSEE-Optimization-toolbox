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
#   N. Gould, private communication.

#   SIF input: N. Gould, Jan 1996

#   classification OUR2-AN-V-0

Translated to scilab from AMPL by Sharvani Laxmi Somayaji and Pranshu Malhotra as a part of FOSSEE internship, 2021
*/

/*-------test case outputs---------
For N=10, in scilab
fval = -9.
 output  = 

  Iterations = 16
  Cpu_Time = 0.048
  Objective_Evaluation = 42
  Dual_Infeasibility = 0.0000003
  Message = "Optimal Solution Found"
For same N, on NEOS server:
objective = -9.0000000000000000e+00 (both scaled and unscaled)

For N=1000, in scilab
fval = -999.
 output  = 

  Iterations = 16
  Cpu_Time = 0.686
  Objective_Evaluation = 34
  Dual_Infeasibility = 0.0000003
  Message = "Optimal Solution Found"
For same N, on NEOS server:
objective = -9.9900000000000000e+02  (both scaled and unscaled)

*/
funcprot(0);
N = 10; 
//N= 1000;
x0 = ones(1,N);
//x0(1,1:N) = 1;

function y=f(x)
    N = 10; 
    //N= 1000;
    y=0;
    for i=(1:(N-1)) 
        y = y + (cos(-0.5.*x(1,i+1) + x(1,i).^2));
    end
endfunction

function g=fGrad(x)
    N = 10;     
    //N= 1000;
    g = zeros(N,1);

    for i=(2:(N-1)) 
        g(i,1) = g(i,1) + (0.5).*(sin(-0.5.*x(i) + x(i-1).^2)) - 2.*x(i).*(sin(-0.5.*x(i+1) + x(i).^2));
    end
    g(1,1) = (sin(0.5.*x(1+1) - x(1).^2)*2).*x(1);
    g(N,1) = (sin(0.5.*x(N) - x(N-1).^2)).*(-0.5);
    /*
    for i=(1:(N-1)) 
        g(i+1,1) = g(i+1,1) + (0.5).*(sin(-0.5.*x(1,i+1) + x(1,i).^2));
    end
    for i=(1:(N-1)) 
        g(i,1) = g(i,1) - 2.*x(1,i).*(sin(-0.5.*x(1,i+1) + x(1,i).^2));
    end
    */
endfunction

options = struct("MaxIter", [10000], "CpuTime", [6000], "GradObj", fGrad, "Hessian","off","GradCon","off","HessianApproximation",[1]);
[x,fval,exitflag,output] =fot_fmincon(f, x0,[],[],[],[],[],[],[],options)
