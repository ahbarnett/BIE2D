function err=checkgrad(f,df,x0)
% CHECKGRAD   verify an function in R^2 has correct analytic gradient
%
% err = checkgrad(f,df) returns error between finite-differencing approx to
%  grad f and the function df. f is a function handle taking a stack of m 3-cpt
%  col vecs and returning a row vec of m values. df only need map a single
%  col vec to a col vec gradient vector (ie not be vectorized).
% err = checkgrad(f,df,x0) enforces the testpoint, which is otherwise random
%  in [-1,1]^2
%
% a self-test is done with no inputs

% Barnett 7/18/16, 2/7/21

if nargin==0, test_checkgrad; return; end
if nargin<3, x0 = rand(2,1)*2 -1 ; end
eps = 1e-5;
x = repmat(x0,[1,5]);  % set up coords to eval f at
for i=1:2, v = zeros(2,1); v(i) = 1;  % unit vec
  x(:,2*i) = x(:,2*i)-eps*v; x(:,2*i+1) = x(:,2*i+1)+eps*v;
end
u = f(x);
du = [u(3)-u(2);u(5)-u(4)]/(2*eps);  % hmm, I guess u(1) ignored
err = norm(du-df(x0));

%%%%%%%%
function test_checkgrad    % throw it a simple pair (f,df) that is correct
a = [1;-1.6];          % the const grad
f = @(x) 7.0 + a(1)*x(1,:) + a(2)*x(2,:);  % linear func
df = @(x) repmat(a,[1 size(x,2)]);
checkgrad(f,df)
checkgrad(f,df,[2;3])
