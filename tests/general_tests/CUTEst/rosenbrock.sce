//--------rosenbrock----------
//Written by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
/*
-------------test case outputs--------------

In scilab
with hessian approximation off:

 fval  = 0.2488974
 output  = 

  Iterations = 104
  Cpu_Time = 0.423
  Objective_Evaluation = 105
  Dual_Infeasibility = 0.0000013
  Message = "Optimal Solution Found"

 x  = 0.5022024   0.2488983
 
with hessian approximation on:

  fval  = 0.2488970
 output  = 

  Iterations = 36
  Cpu_Time = 0.172
  Objective_Evaluation = 42
  Dual_Infeasibility = 0.0000011
  Message = "Optimal Solution Found"

 x  = 0.5022027   0.2488986

*/

funcprot(0);
//initial value
x0=[-1,2];
//linear inequality constraint
A=[1,2];
b=1;
//A =[];
//b=[];
//linear equality constraint
Aeq=[];
beq=[];
//lower and upper bounds
lb=[];
ub=[];

//objective function
function y=f(x)
    y=100*(x(2) - x(1)^2)^2+(1-x(1))^2;
endfunction

//gradient for the objective function
function y = fGrad(x)
    y = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
        200*(x(2)-x(1)^2)];
endfunction

//non linear constraint
function [c,ceq]=nlc(x)
    c=[x(1)^2+x(2)^2-1];
    ceq=[];
endfunction

function y = lHess(x,obj,lambda)
    //Hessian of objective
    H = [1200*x(1)^2-400*x(2)+2, -400*x(1);
                -400*x(1), 200];
    //Hessian of nonlinear inequality constraint
    Hc = [2,0;0,2];
    //lambda.ineqnonlin: The Lagrange multipliers for the nonlinear inequality constraints.
    y = obj*H + lambda(1)*Hc;
    
endfunction

//gradient for the non linear constraint
function [cg, ceqg] = cGrad(x)
    cg = [2*x(1),2*x(2)];
    ceqg = [];
endfunction

//With hessian approximation off
//options = struct("MaxIter", [15000], "CpuTime", [5000], "GradObj", fGrad, "Hessian", lHess,"GradCon", cGrad,"HessianApproximation", [0]);
//With hessian approximation on
options = struct("MaxIter", [15000], "CpuTime", [5000], "GradObj", fGrad, "Hessian", "off","GradCon", cGrad,"HessianApproximation", [1]);

[x,fval,exitflag,output] =fot_fmincon(f, x0,A,b,Aeq,beq,lb,ub,nlc,options)
