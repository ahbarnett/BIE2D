function [x,nsteps] = bisect(f,a,b,tol)
% BISECT  Conduct bisection search to find root of 1D function in interval.
%
% x = bisect(f,a,b) estimates x where f(x) changes sign in [a,b].
%
% x = bisect(f,a,b,tol) controls the absolute tolerance in x to reach
%  (the default is 1e-6).
%
% [x, nsteps] = ... also returns the number of steps used, for diagnosis.
%
% Without arguments, a self-test is done.
%
% Typically O(log_2((b-a)/tol)) function evals are needed

% Barnett 6/19/21, borrowing code from my M56 2014 student exercise.
if nargin==0, test_bisect; return; end
if nargin<4, tol=1e-6; end

maxsteps = 60;
if b<=a, error('b must be greater than a!'); end
if tol<=0, error('tol must be positive!'); end
fa = f(a);
fb = f(b);
if sign(fa)==sign(fb), error('f(a) and f(b) have same sign!'); end
nsteps = 0;
while (b-a)/2 > tol
  x = (a+b)/2;
  nsteps = nsteps+1;
  if nsteps > maxsteps
    warning('too many steps: returning current x');
    break
  end
  fx = f(x);
  if sign(fx)==sign(fa)
    a = x;
    fa = fx;
  else
    b = x;
    fb = fx;
  end
end
x = (a+b)/2;

%%%%%%%%%
function test_bisect
[x,nsteps] = bisect(@sin,3,4,1e-14)
x-pi
%[x,nsteps] = bisect(@(x) x,0,1)
%[x,nsteps] = bisect(@(x) 1-x,0,1)
try, [x,nsteps] = bisect(@sin,3,2); catch me, ['ok: ',me.message], end
try, [x,nsteps] = bisect(@sin,3,4,-1); catch me, ['ok: ',me.message], end
try, [x,nsteps] = bisect(@sin,2,3); catch me, ['ok: ',me.message], end
try, [x,nsteps] = bisect(@(x) 0*x,2,3); catch me, ['ok: ',me.message], end
