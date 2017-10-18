function [f,g,H] = fun(x,fnum)
% Computes the following:
%   f: function scalar value at x
%   J: Jacobian vector (df/dx) at x
%   H: Hessian matrix (d2f/dx2) at x
switch fnum
	case 0
		[f,g,H] = fun0(x);
	case 1
		[f,g,H] = fun1(x);
    case 2
        [f,g,H] = fun2(x);
end 
end

%%
function [f,g,H] = fun0(x)
n = length(x);
g = zeros(n,1);
H = zeros(n,n);
x1 = x(1);
x2 = x(2);
A = 0.4;
a = 0.5;
b = 2.4;
f = A*(a*(x1-1)^2 + b*(x2+2)^2);
g(1) = A*(2*a*(x1-1));
g(2) = A*(2*b*(x2+2));
H(1,1) = A*2*a;
H(1,2) = 0;
H(2,1) = 0;
H(2,2) = A*2*b;
end

%%
function [f,g,H] = fun1(x)
n = length(x);
g = zeros(n,1);
H = zeros(n,n);
x1 = x(1);
x2 = x(2);
A = 0.5;
f = A*(sin(x1) + sin(x2) + 0.1*x1*x1 + 0.1*x1*x2 + 0.1*x2*x2);
g(1) = A*(cos(x1) + 0.2*x1 + 0.1*x2);
g(2) = A*(cos(x2) + 0.1*x1 + 0.2*x2);
H(1,1) = A*(-sin(x1) + 0.2);
H(1,2) = A*(0.1);
H(2,1) = A*(0.1);
H(2,2) = A*(-sin(x2) + 0.2);
end

%%
function  [f,g,H] = fun2(x)
n = length(x);
g = zeros(n,1);
H = zeros(n,n);
x1 = x(1);
x2 = x(2);
f = cos(x1)*cos(x2);
g(1) = - sin(x1)*cos(x2);
g(2) = -sin(x2)*cos(x1);
H(1,1) = -cos(x1)*cos(x2);
H(1,2) = sin(x1)*sin(x2);
H(2,1) = sin(x2)*sin(x1);
H(2,2) = -cos(x2)*cos(x1);
end
