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

#   Source: Nick Gould

#   SIF input: Nick Gould, September 1997.

#   classification SUR2-AN-V-0

Translated to scilab from AMPL by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
*/

/*
---------test case output---------
For N=100,K=1, in scilab
CPU time: 14.156000
fval = -10031.629
For same N, on NEOS server:
fval = -1.0031629024133126e+04

*/

funcprot(0);
N = 100; K = 1;
//N = 10; K = 1;
//N=10000; K =10;

xx = zeros(N,1);
x0 = zeros(N,1);

Q = zeros(N,1);
iterN = 1:N;
x0(iterN) = 0.0001/(N+1);

function y=f(x)
    N=100;K=1;
    //N = 10; K = 1;
    //N=10000; K =10;
    
    Q = zeros(N,1);
    for i=1:N
        if i<=(N-K)
            Q(i) = (sum(x(i:(i+K))));
        else
            Q(i) = (sum(x(i:N)));
        end
    end
    i=1:N;
    y = sum(Q(i).*(Q(i).*(Q(i).^2-20)-0.1));

endfunction

function g=fGrad(x)
    N=100;K=1;
    //N = 10; K = 1;
    //N=10000; K =10;
    
    Q = zeros(N,1);
    for i=1:N
        if i<=(N-K)
            Q(i) = (sum(x(i:i+K)));
        else
            Q(i) = (sum(x(i:N)));
        end
    end
    
    dfdQi = zeros(N,1);
    for i=1:N
        dfdQi(i) = (4.*Q(i).^3)-(40.*Q(i))-0.1;
    end
    g= zeros(N,1);
    for i=1:N
        if i<=(K)
            g(i) = sum(dfdQi(1:i));
        else
            g(i) = sum(dfdQi((i-K):i));
        end
    end
    
endfunction

options = struct("MaxIter", [1000000], "CpuTime", [60000], "GradObj", fGrad, "Hessian","off","GradCon","off","HessianApproximation", [1] );
[x,fval,exitflag,output] =fot_fmincon(f, x0,[],[],[],[],[],[],[],options)

