function [u un uxx uxy uyy] = LapSLP(t,s,dens)
% LAPSLP   Evaluate Laplace single-layer potential from curve to targets
%
% This evaluates the 2D Laplace single-layer potential for the density tau,
%
%   u(x) = (1/2pi) int_gamma log(1/r) tau(y) ds_y,   where r:=x-y,  x,y in R2,
%
%  using the native quadrature rule on the source segment gamma, where point
%  values of tau are given.
%
% [u un] = LapSLP(t,s,dens) evaluates potential and its target-normal
%  derivative.
%
% [u un uxx uxy uyy] = LapSLP(t,s,dens) also returns 2nd derivs of potential.
%
% [A An ...] = LapSLP(t,s) or LapSLP(t,s,[]) returns matrix which maps a
%  density vector to the vector of potentials (A), and possibly target-normal
%  derivatives (An) and Hessians (Axx Axy Ayy).
%
% Tested by: LAPINTDIRBVP
% apparently I should add the sameseg situtation for traction here
% yikes
%
% Crude native quadr and O(NM) RAM for now
% todo: make O(N+M) & incorporate Gary's scf

% Barnett 6/27/16; Jun Wang added 2nd derivs, Oct 2018.
  
if nargout==1
  u = LapSLPmat(t,s);
  if nargin>2 && ~isempty(dens)
    u = u * dens;
  end
elseif nargout==2
  [u un] = LapSLPmat(t,s);
  if nargin>2 && ~isempty(dens)
    u = u * dens;
    un = un * dens;
  end
else
 % just add in a second deriv. for comparison
 % haven't done the self case yet
  [u,un,uxx,uxy,uyy]=LapSLPmat(t,s); 
  if nargin>2 && ~isempty(dens)
    u = u*dens;
    un = un*dens;
    uxx = uxx*dens;
    uxy = uxy*dens;
    uyy = uyy*dens;
  end
end
%%%%%%

function [A An A11 A12 A22] = LapSLPmat(t,s)
% [A An] = LapSLPmat(t,s)
% plain single-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg.
% Kress quadrature used for self-interaction (assumes global quad).
% No jump included on self-interaction of derivative (ie PV integral).

% Barnett 6/12/16 from stuff since 2008; bsxfun speedups 6/28/16
N = numel(s.x);
d = bsxfun(@minus,t.x,s.x.');          % C-# displacements mat
r1= real(d);
r2= imag(d);
ss = sameseg(t,s);
if ss
  A = -log(abs(d)) + circulant(0.5*log(4*sin(pi*(0:N-1)/N).^2));   % peri log
  A(diagind(A)) = -log(s.sp);                       % diagonal limit
  m = 1:N/2-1; Rjn = ifft([0 1./m 2/N 1./m(end:-1:1)])/2;  % Kress Rj(N/2)/4pi
  A = A/N + circulant(Rjn); % includes SLP prefac 1/2pi. Kress peri log matrix L
  A = bsxfun(@times, A, s.sp.');   % do speed factors (2pi/N weights already)
else
  A = bsxfun(@times, log(abs(d)), -(1/2/pi)*s.w(:)');  % prefactor & wei
end
if nargout>=2                      % apply D^T (flips sign and src deriv)
  An = real(bsxfun(@rdivide,(1/2/pi)*(-t.nx),d));    % complex form of dipole
  if ss
    An(diagind(An)) = -s.cur*(1/4/pi);  % self? diagonal term for Laplace
  end
  An = bsxfun(@times, An, s.w(:)');  % quadr wei
end

if nargout>=3
  A11=(1/2/pi)*(r1.^2-r2.^2)./abs(d).^4;
  A12=(1/2/pi)*2*r1.*r2./abs(d).^4;
  A22=(1/2/pi)*(r2.^2-r1.^2)./abs(d).^4;

  A11=bsxfun(@times, A11, s.w(:)');
  A12=bsxfun(@times, A12, s.w(:)');
  A22=bsxfun(@times, A22, s.w(:)');
end
