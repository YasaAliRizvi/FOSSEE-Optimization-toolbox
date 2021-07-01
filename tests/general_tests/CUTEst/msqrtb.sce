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

#   Source:  problem 204 (p. 93) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NQR2-AN-V-V

Translated to scilab from AMPL by Yasa Ali Rizvi as a part of FOSSEE internship, 2021
*/

/*
-------test case outputs----------
//  Bestknown Objective = 0.00000000

for P=3 in scilab,
 fval  = 0.
 output  = 

  Iterations = 90
  Cpu_Time = 0.404
  Objective_Evaluation = 135
  Dual_Infeasibility = 0
  Message = "Optimal Solution Found"

For same P, on NEOS server:
objective = 0.0000000000000000e+00 (both scaled and unscaled)
*/
P = 3;
//P = 4;
N = P^2;
B = zeros(P,P);
A = zeros(P,P);
x0 = zeros(P,P);

for i=1:P
    for j=1:P 
        if (i==3 && j==1) then
            B(i,j) = 0;
        else
            B(i,j) = sin(((i-1)*P + j)^2);
        end
    end
end

for i=1:P
    for j=1:P
        x0(i,j) = 0.2*B(i,j);
    end
end

x0 = matrix(x0',N,1);

function y=f(x)
    y=0;
endfunction

function [c, ceq]=nlc(x)
    P = 3;
    //P = 4;
    N = P^2;

    c = [];
    ceq = [];
    x = (matrix(x,P,P))';

    for i=1:P
        for j=1:P 
            if (i==3 && j==1) then
                B(i,j) = 0;
            else
                B(i,j) = sin(((i-1)*P + j)^2);
            end
        end
    end

    for i=1:P 
        for j=1:P
            sum = 0;
            for k=1:P
                sum = sum + B(i,k)*B(k,j);
            end
            A(i,j) = sum;
        end
    end

    for i=1:P
        for j=1:P
            sum=0;
            for t=1:P
                sum = sum + x(i,t)*x(t,j);
            end
            sum = sum-A(i,j);
            ceq = [ceq, sum];
        end
    end
endfunction

function y=fGrad(x)
    P = 3;
    //P = 4;
    N = P^2;
    y = zeros(1, N);
    //y(1,1:N) = 0;
endfunction

function [cg, ceqg]=cGrad(x)
    P = 3;
    //P = 4;
    N = P^2;
    cg=[];
    nc=1;
    ceqg = zeros(nc, N);
    for i=1:P
        for j=1:P
            for t=1:P
                if (i==t && t==j) then
                    ceqg(nc,(i-1)*P+j) = 2*x((i-1)*P+j);
                else
                    ceqg(nc,(i-1)*P+t) = x((t-1)*P+j);
                    ceqg(nc,(t-1)*P+j) = x((i-1)*P+t);
                end
            end
            nc = nc+1;
        end
    end
endfunction

options = struct("MaxIter", [300000], "CpuTime", [60000],"HessianApproximation", [1], "GradObj", fGrad, "Hessian","off","GradCon", cGrad);
[x, fval, exitflag, output] = fot_fmincon(f, x0,[],[],[],[],[],[],nlc,options)
