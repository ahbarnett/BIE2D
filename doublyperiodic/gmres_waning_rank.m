% test GMRES when increasingly rank-deficient in one eigenvalue.
% Barnett 1/11/18
clear

% build a square A close to rank-deficient
n=1000;
T = randn(n);
[Q,R] = qr(T);  % Q real
lam = 1.0 + rand(n,1);  lam(end) = 1e-12;   %lam(end-1) = 1e-14;
L = diag(lam);                 % symm part, sets eigvals if symm
L = L + (5.0/sqrt(n))*triu(randn(n),1);    % control jordan non-normal, non-symm
%L(1:end-1,end) = 0;   % try reducing last col vals, near-nullspace normal
A = Q'*L*Q;

% explore A
l=eig(A);
S = svd(A); disp('last few sing vals of A:'), S(end-5:end)
fprintf('cond(A)=%.3g,   min and max eig of A: %.3g,  %.3g\n',S(1)/S(end), min(abs(l)),max(abs(l)))
figure; subplot(1,2,1); plot(l,'.'); axis equal; title('eig A')

% set up lin sys
x = randn(n,1);
b = A*x;

% solve it
y=A\b;
fprintf('mldivide: ||resid|| = %.3g, ||y||=%.3g, ||y-x||=%.3g\n',norm(A*y-b),norm(y),norm(y-x))

[z,flag,relres,iter,resvec]=gmres(A,b,[],1e-14,n);
fprintf('gmres: ||resid|| = %.3g, ||z||=%.3g, ||z-x||=%.3g\n',norm(A*z-b),norm(z),norm(z-x))
its=iter(2),
subplot(1,2,2); semilogy(resvec,'+-'); title('GMRES residual norm');

% Conclusions:
% * one small eigval inserts ~20 extra iters, at that resid level.
% * A strict zero eigval has no effect on getting a good resid.


% Some refs, relevant ?
% see: Brown et al 1997, and
%https://scicomp.stackexchange.com/questions/582/stopping-criteria-for-iterative-linear-solvers-applied-to-nearly-singular-system
