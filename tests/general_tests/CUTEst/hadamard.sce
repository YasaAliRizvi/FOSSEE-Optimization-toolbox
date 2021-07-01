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

#   Source:  A suggestion by Alan Edelman (MIT).

#   SIF input: Nick Gould, Nov 1993.

#   classification LQR2-RN-V-V

Translated to scilab from AMPL by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
*/

/*
--------test case outputs----------

For n=4 in scilab, with hessianApproximation off (and passing hessian function)
fval  = 1.0000295
 output  = 

  Iterations = 11
  Cpu_Time = 0.148
  Objective_Evaluation = 14
  Dual_Infeasibility = 0.0000799
  Message = "Optimal Solution Found"
For n=4 in scilab, without passing hessian function and with hessianApproximation off
 fval  = 

   1.0000295
 output  = 

  Iterations = 11
  Cpu_Time = 59.708
  Objective_Evaluation = 14
  Dual_Infeasibility = 0.0000799
  Message = "Optimal Solution Found"
  
  For same n on NEOS server, 
  objective = 1.0000000045464403e+00 (scaled and unscaled both)
  
For n=8 in scilab, without passing hessian function and with hessianApproximation off
 fval  = 

   1.0000058

 output  = 

  Iterations = 12
  Cpu_Time = 1319.887
  Objective_Evaluation = 39
  Dual_Infeasibility = 0.0000091
  Message = "Optimal Solution Found"
For same n on NEOS server, 
  objective = 1.0000001504207570e+00 (scaled and unscaled both)
For same n in gamsworld site list, fval = 1.00000000
  
*/
funcprot(0);

n=4;
//n = 8;

N = 1:n;
maxval = 0;
x0 = zeros(n^2+1,1);
x0(1,1)= maxval;
x0(2:(n^2+1),1) = 1;

//Lower and Upper bounds
lb = zeros(n^2+1);
lb(1) = 0;
lb(2:(n^2+1)) = -1;

ub = zeros(n^2+1);
ub(1) = %inf;
ub(2:(n^2+1)) = 1;

//Linear Inequality constraints
//constraint 1
A1 = zeros(n^2, n^2+1); b1 = zeros(n^2,1);
//constraint 2
A2 = zeros(n^2, n^2+1); b2 = zeros(n^2,1);
for j =2:(n^2+1)
    A1(j-1,j) = 1;
    A2(j-1,j) = -1;
end
i = 1:(n^2);
A1(i,1) = -1;
A2(i,1) = -1;
A = [A1;A2]; b=[b1;b2];

function y=f(x)
     maxval = x(1);
     y = abs(maxval);
endfunction

function g=fGrad(x)
    n=4;
    //n=8;
      g = zeros(1,n^2+1);
     if maxval>=0
          g(1,1) = 1;
     else
          g(1,1) = -1;
     end
endfunction

function [c,ceq]=nlc(x)
    n=4;
    //n=8;
    ceq = zeros(n,n);
    N=1:n;
    k = 1:n;
    //QQT = zeros(n,n);
    //maxval = x(1);
    Q = matrix(x(2:(n^2+1)),n,n);
    for i=N
        for j=N
            //QQT(i,j) = sum(Q(k,i).*Q(k,j));
            //ceq(i, j) = QQT(i,j) - n;
            ceq(i, j) = sum(Q(k,i).*Q(k,j)) - n;
        end
    end
    ceq = matrix(ceq, 1, n^2);
    c=[];
endfunction

function [gc, gceq]=cGrad(x)
    n=4;
    //n=8;
    ceq = zeros(n,n);
    QQT = zeros(n,n);
    N=1:n;
    k = 1:n;
    
    maxval = x(1);
    Q = matrix(x(2:(n^2+1)),n,n);

    gc = [];
    //n^2 constraint equations and n^2+1 variables
    gceq = zeros(n^2, n^2+1);
    Qgrad = zeros(n,n);
    
    for i=1:n
        for j=1:n
                Qgrad = zeros(n,n);
                for k=1:n
                    Qgrad(k, i) = Qgrad(k, i)+ Q(k, j);
                end
                for k =1:n
                    Qgrad(k, j) = Qgrad(k, j)+ Q(k, i);
                end
                k=1:(n^2);
                Qgradtemp = matrix(Qgrad,n^2,1);
                gceq( i+(j-1)*(n), 1) = 0;
                gceq(i+(j-1)*(n), k+1) = (Qgradtemp(k, 1))'; 
        end 
    end
endfunction

function y = lHess(x,obj,lambda)
    n=4;
    //n=8;
    
    //Hessian of objective
    //H = 0;
    //Hessian of nonlinear equality constraint
    //n^2 constraint equations and n^2+1 variables
    Hct = zeros(n^2+1, n^2+1);
    for i=1:n
        for j=1:n
                Qggrad = zeros( n^2+1, n^2+1);
                for k=1:n
                    Qggrad(1+k+(j-1)*n, 1+k+(i-1)*n) = Qggrad(1+k+(j-1)*n, 1+k+(i-1)*n) +1;
                end
                for k=1:n
                    Qggrad(1+k+(i-1)*n, 1+k+(j-1)*n) = Qggrad(1+k+(i-1)*n,1+ k+(j-1)*n) +1; 
                end
                Hct = Hct + lambda(i+(j-1)*n)*Qggrad;
        end 
    end
    y = Hct;
endfunction

//Options
options=struct("MaxIter", [150000], "CpuTime", [50000], "GradObj", fGrad, "Hessian",lHess ,"GradCon", cGrad,"HessianApproximation", [0]);
//Calling Ipopt
[x,fval,exitflag,output] =fot_fmincon(f, x0,A,b,[],[],lb,ub,nlc,options)
