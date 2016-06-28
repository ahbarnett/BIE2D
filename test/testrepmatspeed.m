% repmat is much slower than outer product with ones. Barnett 6/27/16
% Timings on 4-core i7 laptop:

n=300;

tic;a=rand(1,n);for i=1:1e4,b=a*ones(n,1);end,toc
%Elapsed time is 0.041663 seconds.

tic;a=rand(1,n);for i=1:1e4,b=repmat(a,[n 1]);end,toc
%Elapsed time is 1.449605 seconds.

tic;a=rand(1,n);for i=1:1e4,b=bsxfun(@times,a,ones(n,1));end,toc
%Elapsed time is 1.240348 seconds.
