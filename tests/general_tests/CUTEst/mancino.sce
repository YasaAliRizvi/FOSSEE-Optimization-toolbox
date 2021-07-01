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
#   E. Spedicato,
#   "Computational experience with quasi-Newton algorithms for
#   minimization problems of moderate size",
#   Report N-175, CISE, Milano, 1975.

#   See also Buckley #51 (p. 72), Schittkowski #391 (for N = 30)

#   SIF input: Ph. Toint, Dec 1989.
#              correction by Ph. Shott, January, 1995.

#   classification SUR2-AN-V-0

Translated to scilab from AMPL by Yasa Ali Rizvi as a part of FOSSEE internship, 2021
*/

/*
-------test case outputs--------
//Bestknown Objective = 0.00000000; (f = 1.03833e-21;N=40)

//with fminunc
//N = 40; //fval = 2.467D-11; fval(with grad) = 2.830D-11
//N = 20; //fval=   1.668D-12(without grad);  NEOS=1.95804e-15;
//N = 25; //fval=   2.589D-11(without grad);  NEOS=4.30596e-14;
//N = 30; //fval=   2.637D-10(without grad);  NEOS=5.35396e-13;
//N = 35; //fval=   1.005D-10(without grad);  NEOS=4.49868e-12;
//N = 40; //fval=   9.067D-10(without grad);  NEOS=1.03833e-21;
//N = 60; //fval=   1.277D-10(without grad);  NEOS=1.65624e-21;
//N = 100; //fval=   2.406D-09(without grad); NEOS=7.55445e-19;

for N=60 in scilab,
 fval  = 2.300D-09
 output  = 

  Iterations = 9
  Cpu_Time = 4.69
  Objective_Evaluation = 16
  Dual_Infeasibility = 0.1131235
  Message = "Optimal Solution Found"
for same N value on NEOS server,
    objective = 3.2539333553302507e-27 (scaled)
    objective  = 1.6957851316545521e-21(unscaled)
*/

N=60;
//x0 = zeros(N^2+2*N, 1);
x0=zeros(N);

for i=1:N
    sum1 = 0;
    for j=1:N
        sum1 = sum1 + sqrt(i/j)*((sin(log(sqrt(i/j))))^5+(cos(log(sqrt(i/j))))^5);
    end
    x0(i) = -8.710996D-4*((i-50)^3 + sum1);
end

function y=f(x)
// k=1;
// for i=1:N
// for j=1:N
//     v(i,j) = x(N+k);
//     k=k+1;
// end
// end
// alpha = x(N+N^2+1:2*N+N^2);
    N=60;
    v=zeros(N,N);
    alpha=zeros(N,1);
    for i=1:N
        for j=1:N
            v(i, j) = sqrt(x(i)^2 + i/j);
        end
    end

    for i=1:N
        sum2 = 0;
        for j=1:N
            sum2 = sum2 + v(i,j)*((sin(log(v(i,j))))^5 + (cos(log(v(i,j))))^5);
        end
        alpha(i) = 1400*x(i) + (i-50)^3 + sum2;
    end

    sum2 = 0;
    for i=1:N
        sum2 = sum2 + alpha(i)^2;
    end
    y = sum2;
endfunction

function y=fGrad(x)
// k=1;
// for i=1:N
// for j=1:N
//     v(i,j) = x(N+k);
//     k=k+1;
// end
// end
// alpha = x(N+N^2+1:2*N+N^2);
    N=60;
    v = zeros(N,N);
    for i=1:N
        for j=1:N
            v(i, j) = sqrt(x(i)^2 + i/j);
        end
    end

    alpha = zeros(N, 1);
    for i=1:N
        sum1 = 0;
        for j=1:N
            sum1 = sum1 + v(i,j)*((sin(log(v(i,j))))^5 + (cos(log(v(i,j))))^5);
        end
        alpha(i) = 1400*x(i) + (i-50)^3 + sum1;
    end

    dv = zeros(N,N);
    for i=1:N
        for j=1:N
            dv(i,j) = x(i)/sqrt(x(i)^2 + i/j);
        end
    end

    dalpha = zeros(N,1);
    for i=1:N
        sum2 = 0;
        for j=1:N
            sum2 = sum2 + ( dv(i,j)*((sin(log(v(i,j))))^5 + (cos(log(v(i,j))))^5) + v(i,j)*( (5*(sin(log(v(i,j))))^4 * cos(log(v(i,j))) * 1/v(i,j) * dv(i,j)) + (5*(cos(log(v(i,j))))^4 * (-sin(log(v(i,j)))) * 1/v(i,j) * dv(i,j)) ));
        end
        dalpha(i) = 1400 + sum2;
    end

    for i=1:N
        y(1, i) = 2*alpha(i)*dalpha(i);
    end
endfunction

options = struct("MaxIter", [30000], "CpuTime", [6000], "HessianApproximation", [1], "GradObj", fGrad, "Hessian","off","GradCon","off");
[x,fval,exitflag,output] =fot_fmincon(f, x0,[],[],[],[],[],[],[],options)
