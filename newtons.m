function [x,fx,xx] = newtons(f,x0,TolX,MaxIter,varargin)
%newtons.m to solve a set of nonlinear eqs f1(x)=0, f2(x)=0,..
%input: f = 1^st-order vector ftn equivalent to a set of equations
% x0 = the initial guess of the solution
% TolX = the upper limit of |x(k) - x(k - 1)| % MaxIter = the maximum # of iteration
%output: x = the point which the algorithm has reached
% fx = f(x(last))
% xx = the history of x
h = 1e-4; TolFun = eps; EPS = 1e-6;
fx = feval(f,x0,varargin{:});
Nf = length(fx); Nx = length(x0);
if Nf ~= Nx, error('Incompatible dimensions of f and x0!'); end
if nargin < 4, MaxIter = 100; end
if nargin < 3, TolX = EPS; end
xx(1,:) = x0(:).'; %Initialize the solution as the initial row vector
fx0 = norm(fx); %(1)
for k = 1: MaxIter
    dx = -jacob(f,xx(k,:),h,varargin{:})\fx(:);%/;%-[dfdx]?-1*fx
    for l = 1: 3 %damping to avoid divergence %(2)
        dx = dx/2; %(3)
        xx(k + 1,:) = xx(k,:) + dx.';
        fx = feval(f,xx(k + 1,:),varargin{:}); fxn = norm(fx);
        if fxn < fx0, break; end %(4)
    end %(5)
    if fxn < TolFun | norm(dx) < TolX, break; end
    fx0 = fxn; %(6)
end
x = xx(k + 1,:);
if k == MaxIter, fprintf('The best in %d iterations\n',MaxIter), end

function g = jacob(f,x,h,varargin) %Jacobian of f(x)
if nargin < 3, h = 1e-4; end
h2 = 2*h; N = length(x); x = x(:).'; I = eye(N);
for n = 1:N
    g(:,n) = (feval(f,x + I(n,:)*h,varargin{:}) ...
    -feval(f,x - I(n,:)*h,varargin{:}))'/h2;
end
