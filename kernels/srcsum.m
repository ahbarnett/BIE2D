function [A B C] = srcsum(kernel, trlist, phlist, t, s, varargin)
% SRCSUM   Sum any kernel evaluation or matrix over a set of source translations
%
% Note, unlike srcsum2 this handles self-interactions correctly
%
% [A B ...] = srcsum(kernel, trlist, phlist, t, s)
% [A B ...] = srcsum(kernel, trlist, phlist, t, s, param)
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
%  s = wobblycurve(.3,5,100);
%  tau = sin(4*s.t);                            % pick a density
%  u = srcsum(@LapDLP,[-2 0 2], [], s,s, tau);  % do it (incl a self-int)
%  ss = s; v = LapDLP(s,ss,tau);                % check matches direct sum
%  ss.x = s.x-2; v = v + LapDLP(s,ss,tau);
%  ss.x = s.x+2; v = v + LapDLP(s,ss,tau);
%  norm(u-v)
%
% Also see: SRCSUM2

% Barnett 6/12/16, 6/29/16 phlist, self-demo. 6/30/16 "a" field for ext close.

if nargin==0, test_srcsum; return; end

afield = isfield(s,'a');
if isempty(phlist), phlist = 0*trlist+1.0; end
xo = s.x; s.x = xo + trlist(1); if afield, ao = s.a; s.a = ao + trlist(1); end
if nargout==1
  A = kernel(t,s,varargin{:});  % since we don't know the size can't prealloc
  for i=2:numel(trlist)
    s.x = xo + trlist(i); if afield, s.a = ao + trlist(i); end
    A = A + phlist(i)*kernel(t,s,varargin{:});
  end
elseif nargout==2
  [A B] = kernel(t,s,varargin{:});
  for i=2:numel(trlist)
    s.x = xo + trlist(i); if afield, s.a = ao + trlist(i); end
    [Ai Bi] = kernel(t,s,varargin{:});
    A = A + phlist(i)*Ai; B = B + phlist(i)*Bi;
  end
elseif nargout==3
  [A B C] = kernel(t,s,varargin{:});
  for i=2:numel(trlist)
    s.x = xo + trlist(i); if afield, s.a = ao + trlist(i); end
    [Ai Bi Ci] = kernel(t,s,varargin{:});
    A = A + phlist(i)*Ai; B = B + phlist(i)*Bi; C = C + phlist(i)*Ci;
  end
else, error('srcsum cannot handle >3 output args');
end

%%%%%%%%%%%%%%%%
function test_srcsum             % includes self-interactions

% density case...
s = wobblycurve(1,.3,5,100);
tau = sin(4*s.t);                            % pick a density
u = srcsum(@LapDLP,[0 2], [1 3], s,s, tau);  % do it  (include "phase")
ss = s; v = LapDLP(s,ss,tau);                % check matches direct sum
ss.x = s.x+2; v = v + 3*LapDLP(s,ss,tau);   % note phase
norm(u-v)

% matrix case...
u = srcsum(@LapDLP,[0 2], [1 3], s,s);  % do it  (include phase)
ss = s; v = LapDLP(s,ss);                % check matches direct sum
ss.x = s.x+2; v = v + 3*LapDLP(s,ss);   % note phase
norm(u-v)

% 2-cmpt output matrix case...
mu = 0.6;
u = srcsum(@StoDLP,[0 2], [1 3], s,s,mu);  % do it  (include phase)
ss = s; v = StoDLP(s,ss,mu);                % check matches direct sum
ss.x = s.x+2; v = v + 3*StoDLP(s,ss,mu);   % note phase
norm(u-v)

% 2-cmpt multi-output matrix case...
[u p] = srcsum(@StoDLP,[0 2], [1 3], s,s,mu);  % do it  (include phase)
ss = s; [v q] = StoDLP(s,ss,mu);                % check matches direct sum
ss.x = s.x+2; [v2 q2] = StoDLP(s,ss,mu);
v = v + 3*v2; q = q + 3*q2;
norm(u-v)

