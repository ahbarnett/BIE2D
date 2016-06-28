function testrepmatspeed
% repmat timing vs outer product with ones, bsxfun, etc. Barnett 6/27/16
% Timings on 4-core i7 laptop.
% Consider run with:
% profile clear; profile on; testrepmatspeed; profile off; profile viewer

if 1   % basic repeating timing. bsxfun wins over others, slightly
  n=300;
  a=rand(1,n);tic,for i=1:1e4,b=ones(n,1)*a;end,toc
  %Elapsed time is 1.224451 seconds.
  a=rand(1,n);tic, for i=1:1e4,b=repmat(a,[n 1]);end,toc
  %Elapsed time is 1.449605 seconds.
  a=rand(1,n);tic; for i=1:1e4,b=bsxfun(@times,a,ones(n,1));end,toc
  %Elapsed time is 1.240348 seconds.
end

if 1   % mock up bottlenecks in Lap DLP dense mat filling
  n=300; r=1e4; a=rand(n); b=rand(n); d=rand(1,n);
  tic,for i=1:r,c=a./b;end,toc        % is faster! (only 10-20%)
  tic,for i=1:r,c=bsxfun(@rdivide,a,b);end,toc
  tic,for i=1:r,c=repmat(d,[n 1])./a;end,toc   % slow, case where one is rank 1
  %tic,for i=1:r,c=repmat(d',[1 n])./a;end,toc   % case where one is rank 1
  tic,for i=1:r,c=(ones(n,1)*d)./a;end,toc   % beats repmat
  tic,for i=1:r,c=bsxfun(@rdivide,d,a);end,toc   % 2-3x slower - why?
  tic,for i=1:r,c=bsxfun(@times,d,1./a);end,toc   % 4x slower - bad
  tic,b =repmat(d,[n 1]); for i=1:r,c=b./a;end,toc   % as expected, like a./b
  %tic,for i=1:r,c=bsxfun(@times,a,1./d);end,toc  % not always faster
  %tic,for i=1:1e4,c=a.^2.^2;end,toc    % next 2 are same...
  %tic,for i=1:1e4,c=a.*a.*a.*a;end,toc
end
    
if 1  % core of LapDLP timing. Now, for N<1e3, repmat is faster.
  M = 1e3; N = 1e3; reps = 100;
  t.x = rand(M,1); s.x = rand(N,1); s.w = rand(N,1); s.nx = rand(N,1);
  tic, for i=1:reps % repmat version   (8 sec)
  d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
  ny = repmat(s.nx.', [M 1]);      % identical rows given by src normals
  A = (1/2/pi) * real(ny./d);      % complex form of dipole
  %A = A .* repmat(s.w(:)', [M 1]);
  A = bsxfun(@times, A, s.w(:)');
  end, toc
  tic, for i=1:reps % ones version   (14 sec for M = 1e4; N = 1e3; reps = 100)
  d = t.x*ones(1,N) - ones(M,1)*s.x.';   % C-# displacements mat
  ny = ones(M,1)*s.nx.';                 % identical rows given by src normals
  A = (1/2/pi) * real(ny./d);            % complex form of dipole
  A = bsxfun(@times, A, s.w(:)');
  end, toc
  tic, for i=1:reps % bsxfun version  (3.4 sec)
  d = bsxfun(@minus,t.x,s.x.');          % C-# displacements mat
  A = real(bsxfun(@rdivide,(1/2/pi)*s.nx.',d));     % complex form of dipole
  A = bsxfun(@times, A, s.w(:)');
  end, toc
end
% concl: for smallish N, repmat is faster than ones. For col-scaling, bsxfun
%  is faster than repmat. So, leave repmats alone, apart from col-scaling
%  and division change to bsxfun.
