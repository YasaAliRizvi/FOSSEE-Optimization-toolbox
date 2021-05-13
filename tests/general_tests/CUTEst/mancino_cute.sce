//	---Question----

//	Bestknown Objective = 0.00000000; (f = 1.03833e-21;N=40)

//	param N := 100;

//	var x{i in 1..N} := -8.710996D-4*((i-50)^3 + 
//			    sum {j in 1..N} sqrt(i/j)*((sin(log(sqrt(i/j))))^5+ 
//			    (cos(log(sqrt(i/j))))^5));
//	var v{i in 1..N, j in 1..N} = sqrt (x[i]^2 +i/j);
//	var alpha{i in 1..N} = 1400*x[i] + (i-50)^3 + 
//			       sum {j in 1..N} v[i,j]*((sin(log(v[i,j])))^5 + 
//			       (cos(log(v[i,j])))^5);

//	minimize f:
//	sum {i in 1..N} alpha[i]^2;

//	solve;
//	display f;
//	display x;

//with fminunc
//N = 40; //fval = 2.467D-11; fval(with grad) = 2.830D-11
//N = 20; //fval=   1.668D-12(without grad);  NEOS=1.95804e-15;
//N = 25; //fval=   2.589D-11(without grad);  NEOS=4.30596e-14;
//N = 30; //fval=   2.637D-10(without grad);  NEOS=5.35396e-13;
//N = 35; //fval=   1.005D-10(without grad);  NEOS=4.49868e-12;
//N = 40; //fval=   9.067D-10(without grad);  NEOS=1.03833e-21;
//N = 60; //fval=   1.277D-10(without grad);  NEOS=1.65624e-21;
//N = 100; //fval=   2.406D-09(without grad); NEOS=7.55445e-19;

N=100;
//x0 = zeros(N^2+2*N, 1);
x0=zeros(N);

for i=1:N
	sum = 0;
	for j=1:N
		 sum = sum + sqrt(i/j)*((sin(log(sqrt(i/j))))^5+(cos(log(sqrt(i/j))))^5);
	end
	x0(i) = -8.710996D-4*((i-50)^3 + sum);
end

function y=f(x)
	// k=1;
	// for i=1:N
	// 	for j=1:N
	// 		v(i,j) = x(N+k);
	// 		k=k+1;
	// 	end
	// end
	// alpha = x(N+N^2+1:2*N+N^2);
	
	v=zeros(N,N);
	alpha=zeros(N,1);
	for i=1:N
		for j=1:N
			v(i, j) = sqrt(x(i)^2 + i/j);
		end
	end

	for i=1:N
		sum = 0;
		for j=1:N
			sum = sum + v(i,j)*((sin(log(v(i,j))))^5 + (cos(log(v(i,j))))^5);
		end
		alpha(i) = 1400*x(i) + (i-50)^3 + sum;
	end

	sum = 0;
	for i=1:N
		sum = sum + alpha(i)^2;
	end
	y = sum;
endfunction

function y=fGrad(x)
	// k=1;
	// for i=1:N
	// 	for j=1:N
	// 		v(i,j) = x(N+k);
	// 		k=k+1;
	// 	end
	// end
	// alpha = x(N+N^2+1:2*N+N^2);
	v = zeros(N,N);
	for i=1:N
		for j=1:N
			v(i, j) = sqrt(x(i)^2 + i/j);
		end
	end

	alpha = zeros(N, 1);
	for i=1:N
		sum = 0;
		for j=1:N
			sum = sum + v(i,j)*((sin(log(v(i,j))))^5 + (cos(log(v(i,j))))^5);
		end
		alpha(i) = 1400*x(i) + (i-50)^3 + sum;
	end

	dv = zeros(N,N);
	for i=1:N
		for j=1:N
			dv(i,j) = x(i)/sqrt(x(i)^2 + i/j);
		end
	end

	dalpha = zeros(N,1);
	for i=1:N
		sum = 0;
		for j=1:N
			sum = sum + ( dv(i,j)*((sin(log(v(i,j))))^5 + (cos(log(v(i,j))))^5) + ..
			v(i,j)*( (5*(sin(log(v(i,j))))^4 * cos(log(v(i,j))) * 1/v(i,j) * dv(i,j)) + ..
			(5*(cos(log(v(i,j))))^4 * (-sin(log(v(i,j)))) * 1/v(i,j) * dv(i,j)) ));
		end
		dalpha(i) = 1400 + sum;
	end

	for i=1:N
		y(1, i) = 2*alpha(i)*dalpha(i);
	end
endfunction

// options = list("MaxIter", [30000], "CpuTime", [6000], "HessianApproximation", [1], "GradObj", fGrad, "Hessian","off");
//[x,fval, exitflag] = fot_fminunc(f, x0)

options = struct("MaxIter", [30000], "CpuTime", [6000], "HessianApproximation", [1], "GradObj", fGrad, "Hessian","off","GradCon","off");
[x,fval, exitflag] = fot_fmincon(f, x0, [],[],options)