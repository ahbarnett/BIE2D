function [u un] = LapDLP(t,s,dens)
% LAPDLP   Evaluate Laplace double-layer potential from curve to targets
%
% This evaluates the 2D Laplace double-layer potential for the density tau,
%
%   u(x) = (1/2pi) int_gamma (n_y.(x-y))/r^2 tau(y) ds_y,
%   where r:=x-y,  x,y in R2,
%
%  using the native quadrature rule on the source segment, where point
%  values of tau are given.
%
% [u un] = LapDLP(t,s,dens) evaluates potential and its target-normal
%  derviative.
%
% [A An] = LapDLP(t,s) or LapDLP(t,s,[]) returns matrix which maps a
%  density vector to the vector of potentials (A) and target-normal derivatives
%  (An).
%
% Tested by: LAPINTDIRBVP

% Crude native quadr and O(NM) RAM for now. Obviously C/Fortran would not
%  form matrices for the density eval case, rather direct sum w/o/ wasting RAM.
% todo: make O(N+M) & incorporate Gary's scf

% Barnett 6/12/16. Interface change 6/27/16. 
  
if nargout==1
  u = LapDLPmat(t,s);            % local matrix filler
  if nargin>2 && ~isempty(dens)
    u = u * dens;
  end
else
  [u un] = LapDLPmat(t,s);
  if nargin>2 && ~isempty(dens)
    u = u * dens;
    un = un * dens;
  end
end
%%%%%%

function [A An] = LapDLPmat(t,s)
% [A An] = LapDLPmat(t,s)
% plain double-layer kernel matrix & targ n-deriv
% t = target seg (x,nx cols), s = src seg.
% No jump included on self-interaction (ie principal-value integral).
% Self-evaluation for the hypersingular An currently gives inf.

% Barnett 6/12/16 from stuff since 2008. speed repmat->bsxfun 6/28/16
N = numel(s.x);
d = bsxfun(@minus,t.x,s.x.');                 % C-# displacements mat
% mult by identical rows given by source normals...
A = real(bsxfun(@rdivide,(1/2/pi)*s.nx.',d)); % complex form of dipole
if sameseg(t,s)                  % self?
  A(diagind(A)) = -s.cur*(1/4/pi); end        % diagonal term for Laplace
A = bsxfun(@times, A, s.w(:)');
if nargout>1     % deriv of double-layer. Not correct for self-interaction.
  csry = bsxfun(@times, conj(s.nx.'), d);     % (cos phi + i sin phi).r
  % identical cols given by target normals...
  csrx = bsxfun(@times, conj(t.nx), d);       % (cos th + i sin th).r
  r = abs(d);                             % dist matrix R^{MxN}
  An = -real(csry.*csrx)./((r.^2).^2);    % divide is faster than bxsfun here
  An = bsxfun(@times, An, (1/2/pi)*s.w(:)');   % prefac & quadr wei
end
