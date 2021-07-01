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

#   classification OLR2-AN-V-V

Translated to scilab from AMPL by Sharvani Laxmi Somayaji as a part of FOSSEE internship, 2021
*/

/*
--------test case output---------
For n = 20; nx = 2; ny =4; in scilab:
CPU time: 0.112
fval = 0.1126860
For same values, on NEOS server:
fval = 1.1268596668403760e-01

For n= 6;nx=3;ny=3; 
on AMPL(minos solver)
f= 0.1554865896
*/

funcprot(0);

//n = 1000;nx = 5;ny = 10;
n = 20; nx = 2; ny =4;
b = zeros(ny,nx);
for i=1:ny
    for j=1:nx
        b(i,j) = (i-j)/(nx+ny);
    end
end

x0 = zeros((n-1)*nx+n*ny,1);

//linear inequality constraints
//constraint 1
A1 = zeros(n-1, (n-1)*nx+n*ny);
b1=zeros(n-1,1);
i=1:nx;
for t=1:(n-1)
    A1xx = zeros(n-1,nx);
    A1xy = zeros(n,ny);
    A1xy(t,1)=0.5;
    A1xy(t,2)=0.25;
    A1xy(t+1,1)=-1;
    A1xx(t,i) = b(1,i);
    A1xx = matrix(A1xx,1, (n-1)*nx);
    A1xy = matrix(A1xy,1, (n)*ny);
    A1(t,1:((n-1)*nx+n*ny) ) = [A1xx,A1xy];
end
//constraint 2
A2 = zeros((n-1)*(ny-2), (n-1)*nx+n*ny);
b2=zeros((n-1)*(ny-2),1);
i=1:nx;
for t=1:(n-1)
    for j=2:(ny-1)
        A2xx = zeros(n-1,nx);
        A2xy = zeros(n,ny);
        A2xy(t+1,j)=-1;
        A2xy(t,j)=0.5;
        A2xy(t,j-1)=-0.25;
        A2xy(t,j+1)=0.25;
        A2xx(t,i) = b(j,i);
        A2xx = matrix(A2xx,1, (n-1)*nx);
        A2xy = matrix(A2xy,1, (n)*ny);
        A2(t+(j-2)*(n-1),1:((n-1)*nx+n*ny)) = [A2xx,A2xy];
    end
end
  
//constraint 3
A3 = zeros(n-1, (n-1)*nx+n*ny);
b3=zeros(n-1,1);
i=1:nx;
for t=1:(n-1)
    A3xx = zeros(n-1,nx);
    A3xy = zeros(n,ny);
    A3xy(t,ny)=0.5;
    A3xy(t,ny-1)=-0.25;
    A3xy(t+1,ny)=-1;
    A3xx(t,i) = b(ny,i);
    A3xx = matrix(A3xx,1, (n-1)*nx);
    A3xy = matrix(A3xy,1, (n)*ny);
    A3(t,1:((n-1)*nx+n*ny) ) = [A3xx,A3xy];
end

//constraint 4
//Fixing a value:
A4 = zeros(ny, (n-1)*nx+n*ny);
b4=zeros(ny,1);
for i=1:(ny)
    A4xx = zeros(n-1,nx);
    A4xy = zeros(n,ny);
    A4xy(1,i)=1;
    A4xx = matrix(A4xx,1, (n-1)*nx);
    A4xy = matrix(A4xy,1, (n)*ny);
    A4(i,1:((n-1)*nx+n*ny)) = [A4xx,A4xy];
end

Aeq = [A1;A2;A3;A4];
beq = [b1;b2;b3;b4];

function y=f(x)
    //n = 1000;nx = 5;ny = 10;
    n = 20; nx = 2; ny =4;
    //n= 6;nx=3;ny=3;
    i = 1:((n-1)*nx);
    j = ((n-1)*nx+1):((n-1)*nx+n*ny);
    y= sum(sum((x(i)+0.5).^4)) + sum(sum((x(j)+0.25).^4));

endfunction

function g=fGrad(x)
    //n = 1000;nx = 5;ny = 10;
    n = 20; nx = 2; ny =4;
    //n= 6;nx=3;ny=3;
    i = 1:((n-1)*nx);
    j = ((n-1)*nx+1):((n-1)*nx+n*ny);
    g = zeros((n-1)*nx+n*ny , 1);
    g(i,1) = 4.*((x(i)+0.5).^3);
    g(j,1) = 4.*((x(j)+0.25).^3);

endfunction

options = struct("MaxIter", [1000000], "CpuTime", [60000], "GradObj", fGrad, "Hessian","off","GradCon","off","HessianApproximation", [1] );
[x,fval,exitflag,output] =fot_fmincon(f,x0,[],[],Aeq,beq,[],[],[],options)
