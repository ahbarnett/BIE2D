function [A B C] = srcsum2(kernel, trlist, phlist, t, s, varargin)
% SRCSUM2   Sum a kernel eval or matrix over source translations via single call
%
% This is a variant of srcsum that sums over targets in a single kernel call,
%  instead of summing over sources with multiple calls. This is useful for
%  periodized close evaluation. NOTE: it cannot be used for self-interactions!
%
% [A B ...] = srcsum2(kernel, trlist, phlist, t, s)
% [A B ...] = srcsum2(kernel, trlist, phlist, t, s, param)
%
% Inputs:
%  kernel : function of form A = kernel(t,s,...) or [A B] = kernel(t,s,...)
%  trlist : list of complex numbers giving source translations
%  phlist : list of complex numbers giving source phase factors (1 if empty)
%  t      : target segment struct
%  s      : source segment struct
% Outputs:
%  A (& optionally B, C): outputs as from kernel but summed over the source
%
% Example usage:
%  s = wobblycurve(.3,5,100); t = wobblycurve(.3,5,50); t.x = t.x/2; % t neq s
%  tau = sin(4*s.t);                             % pick a density
%  u = srcsum2(@LapDLP,[-2 0 2], [], s,s, tau);  % do it
%  ss = s; v = LapDLP(s,ss,tau);                 % check matches direct sum
%  ss.x = s.x-2; v = v + LapDLP(s,ss,tau);
%  ss.x = s.x+2; v = v + LapDLP(s,ss,tau);
%  norm(u-v)
%
% Also see: SRCSUM

% Barnett 6/30/16 when realized close eval prefers more targs but single src.
% 9/28/17 fixed tt.nx omission.

if nargin==0, test_srcsum2; return; end

if isempty(phlist), phlist = 0*trlist+1.0; end
M = numel(t.x); n = numel(trlist);
tt.x = [];                                 % will store all targets with copies
for i=1:n, tt.x = [tt.x; t.x - trlist(i)]; end  % all transl trgs, by invariance
if isfield(t,'nx'), tt.nx = repmat(t.nx,[n 1]); end  % all trg normals
% (we assumed x and nx are only fields relevant for targets)
if nargout==1
  At = kernel(tt,s,varargin{:});
  A = sumblk(At,M,phlist);
elseif nargout==2
  [At Bt] = kernel(tt,s,varargin{:});
  A = sumblk(At,M,phlist); B = sumblk(Bt,M,phlist);
elseif nargout==3
  [At Bt Ct] = kernel(tt,s,varargin{:});
  A = sumblk(At,M,phlist); B = sumblk(Bt,M,phlist); C = sumblk(Ct,M,phlist);
else
  error('srcsum2 cannot handle >3 output args');
end

function A = sumblk(At,M,phlist)   % crush At along target image sum
n = numel(phlist);
% d will be how many vector components per targ pt (presumably 1 or 2)...
d = size(At,1)/M/n; if d~=1 && d~=2, error('sumblk cant handle # rows At'); end
A = zeros(M*d,size(At,2));
if d==1, ii = 1:M; else, ii = [1:M,n*M+(1:M)]; end    % only handles d=1,2
for i=1:n, A = A + phlist(i)*At(ii+(i-1)*M,:); end


%%%%%%%%%%%%%%%%
function test_srcsum2             % tests non-self interactions only

% density case...
s = wobblycurve(1,.3,5,100); t = wobblycurve(1,.3,5,50); t.x = t.x/2; % t neq s
tau = sin(4*s.t);                            % pick a density
u = srcsum2(@LapDLP,[0 2], [1 3], t,s, tau);  % do it  (include "phase")
ss = s; v = LapDLP(t,ss,tau);                % check matches direct sum
ss.x = s.x+2; v = v + 3*LapDLP(t,ss,tau);   % note phase
norm(u-v)

% matrix case...
u = srcsum2(@LapDLP,[0 2], [1 3], t,s);  % do it  (include phase)
ss = s; v = LapDLP(t,ss);                % check matches direct sum
ss.x = s.x+2; v = v + 3*LapDLP(t,ss);   % note phase
norm(u-v)

% 2-cmpt output matrix case...
mu = 0.6;
u = srcsum2(@StoDLP,[0 2], [1 3], t,s,mu);  % do it  (include phase)
ss = s; v = StoDLP(t,ss,mu);                % check matches direct sum
ss.x = s.x+2; v = v + 3*StoDLP(t,ss,mu);   % note phase
norm(u-v)

% 2-cmpt multi-output matrix case...
[u p] = srcsum2(@StoDLP,[0 2], [1 3],t,s,mu);  % do it  (include phase)
ss = s; [v q] = StoDLP(t,ss,mu);                % check matches direct sum
ss.x = s.x+2; [v2 q2] = StoDLP(t,ss,mu);
v = v + 3*v2; q = q + 3*q2;
norm(u-v)

