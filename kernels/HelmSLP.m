function [u un] = HelmSLP(k,t,s,dens)
% HELMSLP   Evaluate Helmholtz single-layer potential from curve to targets
%
% This evaluates the 2D Helmholtz single-layer potential for the density tau,
%
%   u(x) = (i/4) int_gamma H^{(1)}_0(kr) tau(y) ds_y,  where r:=x-y,  x,y in R2,
%
%  using the native quadrature rule on the source segment gamma, where point
%  values of tau are given. k is the wavenumber (>0).
%
% [u un] = HelmSLP(k,t,s,dens) evaluates potential and its target-normal
%  derivative.
%
% [A An ...] = HelmSLP(k,t,s) or HelmSLP(k,t,s,[]) returns matrix which maps a
%  density vector to the vector of potentials (A), and possibly target-normal
%  derivatives (An)
%
% If t is the same segment as s, the Kress rule for self-evaluation is used,
%  which assumes s is a global periodic trapezoid rule.
%
% Tested by: calling without arguments, and by TESTGRFLAPHELM

% Crude O(NM) RAM for now. Barnett 2/7/21.
if nargin==0, test_HelmSLP; return; end

if nargout==1
  u = HelmSLPmat(k,t,s);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
  end
elseif nargout==2
  [u un] = HelmSLPmat(k,t,s);
  if nargin>3 && ~isempty(dens)
    u = u * dens;
    un = un * dens;
  end
end
%%%%%%

function [A An] = HelmSLPmat(k,t,s)
% [A An] = HelmSLPmat(k,t,s)
% plain single-layer kernel matrix & targ n-deriv, for Helmholtz in R2.
% k = positive wavenumber, t = target seg (x,nx cols), s = src seg.
% Kress quadrature used for self-interaction (assumes global quad) in A, An.
% No jump included on self-interaction of derivative (ie PV integral).
N = numel(s.x);
d = bsxfun(@minus,t.x,s.x.');          % C-# displacements mat
r = abs(d);
ss = sameseg(t,s);
if ~ss                                 % plain far-field rule
  A = bsxfun(@times, besselh(0,k*r), (1i/4)*s.w(:)');  % prefactor & wei
else                                  % Martensen-Kussmaul from MPSpack
  S1 = triu(besselj(0,k*triu(r,1)),1);   % use symm of r matrix for speed
  S1 = -(1/4/pi)*(S1.'+S1);     % next fix it as if diag(r) were 0
  S1(diagind(S1)) = -(1/4/pi);  % S1=M_1/2 of Kress w/o speed fac
  A = triu(besselh(0,k*triu(r,1)),1);         % use symm for speed (arg=0 fast)
  A = (1i/4)*(A.'+A) - S1.*circulant(log(4*sin(pi*(0:N-1)/N).^2)); % A=D2=M_2/2 w/o speed fac
  eulergamma = -psi(1);         % now set diag vals Kress M_2(t,t)/2
  A(diagind(A)) = 1i/4 - eulergamma/2/pi - log((k*s.sp).^2/4)/4/pi;
  %figure; imagesc(real(A)); colorbar; % diag matches (smooth)?
  m = 1:N/2-1; Rjn = -2*pi*ifft([0 1./m 2/N 1./m(end:-1:1)]);  % Kress Rj(N/2)
  A = bsxfun(@times, circulant(Rjn).*S1 + (2*pi/N)*A, s.sp(:)');
end
if nargout>=2                      % apply D^T (flips sign and src deriv)
  if ~ss
    cosker = -real(bsxfun(@times,conj(t.nx),d) ./ r); % -cos ang targ nor to d
    An = besselh(1,k*r) .* cosker;
  else
    % *** todo: D^T with Kress rule?
  end
  An = bsxfun(@times, An, (1i*k/4)*s.w(:)');  % prefactor & wei
end


%%%%
function test_HelmSLP         % only tests pot (ie, value), against reference
b = wobblycurve(1,0.3,5,400); % curve struct
k = 10;                       % wavenumber
t.x = 1.5+1i;                 % far
densfun = @(t) 1+sin(2+4*t+cos(t)) + 1i*sin(1+3*t);  % must be 2pi-per
ker = @(varargin) HelmSLPpotker(k,varargin{:});      % hide the 1st arg (k)
ua = lpevaladapt(t.x,ker,densfun,b,1e-12);
dens = densfun(b.t);          % eval dens at nodes
u = HelmSLP(k,t,b,dens);
disp(max(abs(u-ua)))
