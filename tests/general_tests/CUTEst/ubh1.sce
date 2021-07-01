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

#   Source: unscaled problem 1 
#   (ODE = 1, CLS = 1, GRD = 1, MET = T, SEED = 0.) in
#   J.T. Betts and W.P. Huffman,
#   "Sparse Nonlinear Programming Test Problems (Release 1.0)",
#   Boeing Computer services, Seattle, July 1993.

#   SIF input: Ph.L. Toint, October 1993.

#   classification QLR2-MN-V-V

Translated to scilab from AMPL by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
*/

/*
----------test case outputs-----------

for n=5;t0=0;tf=1000;
    with hessian approximation turned on, no fval was found
 output  = 

  Iterations = 734
  Cpu_Time = 8.157
  Message = "Restoration Failed!"
    
    with hessian approximation turned off and hessian is not provided,  fval  = 1.2444444

 output  = 

  Iterations = 5
  Cpu_Time = 71.808
  Objective_Evaluation = 6
  Dual_Infeasibility = 9.207D-14
  Message = "Optimal Solution Found"
  
  When Hessian approximation is turned off and hessian is provided:
   fval  = 1.2444444
 output  = 

  Iterations = 19
  Cpu_Time = 0.095
  Objective_Evaluation = 20
  Dual_Infeasibility = 0.0000006
  Message = "Optimal Solution Found"
  
for same n and other values, the objective given by ampl(minos): 1.244444444
    
    result obtained from NEOS server for same n and other parameters:
        objective = 1.2444444444444447e+00 (scaled and unscaled both)
    
    
    
for n=2000, t0 =0, tf=1000
    from NEOS server
objecive =    1.1160008156948149e+00 (both scaled and unscaled)
however on gamsworld site, bestknown objective is 14.98048932
*/
funcprot(0);


//n = 2000;
n =5;
tf=1000;

t0 = 0.0;
k = (tf-t0)/n;
x0 = zeros(1, 9*(n+1));
    
//Upper and lower bounds
ub = %inf + zeros(1, 9*(n+1));
lb = -%inf + zeros(1, 9*(n+1));
ub(1, (6*(n+1)+1):9*(n+1)) = 1;
lb(1, (6*(n+1)+1):9*(n+1)) = -1;

//constraint 1
Aeq1 = zeros(3*n, 9*(n+1));
beq1 = zeros(3*n,1);
for t=2:(n+1)
    for i = 1:3
        xytemp =zeros(6,n+1);
        xutemp = zeros(3, n+1);
        xytemp(i,t)=1; xytemp(i,t-1) = -1;
        xytemp(i+3,t-1) = -k/2; xytemp(i+3,t) = -k/2;
        xytemp = matrix(xytemp,1,6*(n+1));
        xutemp = matrix(xutemp,1,3*(n+1));
        Aeq1((t-2)*3+i, 1:(9*(n+1))) = [xytemp,xutemp];
    end
end

//constraint 2
Aeq2 = zeros(3*n, 9*(n+1));
beq2 = zeros(3*n,1);
for t=2:(n+1)
    for i = 1:3
        xytemp =zeros(6,n+1); xutemp = zeros(3, n+1);
        xytemp(i+3,t)=1; xytemp(i+3,t-1) = -1;
        xutemp(i,t-1) = -k/2; xutemp(i,t) = -k/2;
        xytemp = matrix(xytemp,1,6*(n+1));
        xutemp = matrix(xutemp,1,3*(n+1));
        Aeq2((t-2)*3+i, 1:(9*(n+1))) = [xytemp,xutemp];
    end
end

//constraint 3
//Fixing variables
Aeq3 = zeros(12, 9*(n+1));
beq3 = zeros(12,1);
beq3(1,1) = 1000.0;
beq3(2,1) = 1000.0;
beq3(3,1) = 1000.0;
beq3(4,1) = -10.0;
beq3(5,1) = 10.0;
beq3(6,1) = -10.0;
for i = 1:6
    Aeq3(i,i) =1;
end
for i = 1:6
    Aeq3(12-i+1,6*(n+1)-i+1) =1;
end

Aeq = [Aeq1;Aeq2;Aeq3]; beq = [beq1;beq2;beq3];

function y=f(x)
    //n = 2000;
    n =5;
    t0 = 0.0;
    tf=1000;
    k = (tf-t0)/n;
    xu = matrix(x( (6*(n+1)+1): (9*(n+1)) ), 3, n+1);
    
    y=0;
    t = 2:(n+1);
    for i=1:3
        y = y + sum(sum( (k/2).* (xu(i,t).^2 + xu(i,t-1).^2) ) );
    end
endfunction

function g=fGrad(x)
    //n = 2000;
    n =5;
    t0 = 0.0;
    tf=1000;
    k = (tf-t0)/n;
    xu =x((6*(n+1)+1): (9*(n+1)) );
    xu = matrix( xu, 3, n+1);
    gxu = zeros(3, n+1);
    t=2:n;
    for i=1:3
        gxu(i,t) = 2*k.*( xu(i,t) );
        gxu(i,1) =  k.*( xu(i,1) );
        gxu(i,n+1) = k.*( xu(i,n+1) );
    end
    gxu = matrix(gxu,1, 3*(n+1) );
    gxy = zeros(1, 6*(n+1));
    g = [gxy,gxu];
endfunction

function y = lHess(x,obj,lambda)
    n =5;
    t0 = 0.0;
    tf=1000;
    k = (tf-t0)/n;
    //Hessian of objective
    H = zeros(9*(n+1),9*(n+1));
    for i = (6*(n+1)+1): (9*(n+1))
        H(i,i) = 2*k;
    end
    H(6*(n+1),6*(n+1)) = k; H(9*(n+1),9*(n+1)) = k;
    y = obj*H;
endfunction

options = struct("MaxIter", [10000], "CpuTime", [6000], "GradObj",fGrad, "Hessian",lHess,"GradCon","off","HessianApproximation",[0]);
[x,fval,exitflag,output] =fot_fmincon(f, x0,[],[],Aeq,beq,lb,ub,[],options)
