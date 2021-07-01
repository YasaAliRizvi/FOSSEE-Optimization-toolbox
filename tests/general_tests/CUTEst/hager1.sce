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

#   Source: problem P1 in
#   W.W. Hager,
#   "Multiplier Methods for Nonlinear Optimal Control",
#   SIAM J. on Numerical Analysis 27(4): 1061-1080, 1990.

#   SIF input: Ph. Toint, March 1991.

#   classification SLR2-AN-V-V

Translated to scilab from AMPL by Yasa Ali Rizvi as a part of FOSSEE internship, 2021
*/

/*
------test case outputs---------

//Bestknown Objective = 1.31145470; =0.880797(for N=1000)

//N=500; //fval=   0.8807972
//N=700; //fval=   0.8807971 

For N=100, in scilab

 fval  = 0.8807988
 output  = 

  Iterations = 6
  Cpu_Time = 0.132
  Objective_Evaluation = 7
  Dual_Infeasibility = 0.0000002
  Message = "Optimal Solution Found"
For same N, on NEOS server:
fval = 8.8079882788611530e-01 (both scaled and unscaled)

For N=1000, in scilab

 fval  = 0.8807971
 output  = 

  Iterations = 5
  Cpu_Time = 79.41
  Objective_Evaluation = 6
  Dual_Infeasibility = 0.0000004
  Message = "Optimal Solution Found"

For same N, on NEOS server:
fval = 8.8079709547681406e-01 (both scaled and unscaled)

*/
//N = 1000;
N=100;
x0 = zeros(2*N+1, 1);

//Linear equality constraints
beq = zeros(N+1, 1);
beq(1) = 1;
//Aeq = [1];
Aeq = ones(N+1,2*N+1);

Aeq(1, 2:2*N+1) = 0;

for i=2:N+1
	if i>=3 then
		Aeq(i, 1:i-2) = 0;
		Aeq(i, N+2:N+i-1) = 0;
	end
	if i>=2 then
		Aeq(i, i-1) = -N-0.5;
	end
	Aeq(i, i) = N-0.5;
	Aeq(i, i+1:N+1) = 0;
	Aeq(i, N+i) = -1;
	Aeq(i, N+i+1:2*N+1) = 0;
end

function y=f(x)
    //N = 1000;
    N=100;
	sumt = 0 ;
	for i=1:N 
		sumt = sumt + ((x(N+1+i)^2)/(2*N));
	end
	y = 0.5*x(N+1)^2 + sumt;
endfunction

function y=fGrad(x)
    //N = 1000;
    N=100;
    y = zeros(1,2*N+1);
	y(1,1:N) = 0;
	y(1,N+1) = x(N+1);
	for i=1:N
		y(1,N+1+i) = x(N+1+i)/N;
	end
endfunction

options = struct("MaxIter", [3000], "CpuTime", [600], "HessianApproximation", [1], "GradObj", fGrad, "Hessian","off","GradCon","off" );
[x,fval, exitflag, output] = fot_fmincon(f, x0, [], [], Aeq, beq, options)
