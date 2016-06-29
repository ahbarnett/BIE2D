function [A B] = srcsum(kernel, trlist, phlist, t, s, varargin)
% SRCSUM   Sum any kernel evaluation matrix over a set of source translations
%
% [A B] = srcsum(kernel, trlist, phlist, t, s)
% [A B] = srcsum(kernel, trlist, phlist, t, s, param)
%
% Inputs:
%  kernel : function of form A = kernel(t,s,...) or [A B] = kernel(t,s,...)
%  trlist : list of complex numbers giving source translations
%  phlist : list of complex numbers giving source phase factors (1 if empty)
%  t      : target segment struct
%  s      : source segment struct
% Outputs:
%  A (& optionally B): outputs as from kernel but summed over the source

% Example usage:
%  s = wobblycurve(.3,5,100);
%  tau = sin(4*s.t);                            % pick a density
%  u = srcsum(@LapDLP,[-2 0 2], [], s,s, tau);  % do it
%  ss = s; v = LapDLP(s,ss,tau);                % check matches direct sum
%  ss.x = s.x-2; v = v + LapDLP(s,ss,tau);
%  ss.x = s.x+2; v = v + LapDLP(s,ss,tau);
%  norm(u-v)

% Barnett 6/12/16, 6/29/16 phlist, self-demo
% todo: * 3 outputs case

if nargin==0, demo_srcsum; return; end

if isempty(phlist), phlist = 0*trlist+1; end
xo = s.x; s.x = xo + trlist(1);
if nargout==1
  A = kernel(t,s,varargin{:});  % since we don't know the size can't prealloc
  for i=2:numel(trlist)
    s.x = xo + trlist(i);
    A = A + phlist(i)*kernel(t,s,varargin{:});
  end
else
  [A B] = kernel(t,s,varargin{:});
  for i=2:numel(trlist)
    s.x = xo + trlist(i);
    [Ai Bi] = kernel(t,s,varargin{:});
    A = A + phlist(i)*Ai; B = B + phlist(i)*Bi;
  end
end

function demo_srcsum
s = wobblycurve(.3,5,100);
tau = sin(4*s.t);                            % pick a density
u = srcsum(@LapDLP,[-2 0 2], [1 1 1i], s,s, tau);  % do it
ss = s; v = LapDLP(s,ss,tau);                % check matches direct sum
ss.x = s.x-2; v = v + LapDLP(s,ss,tau);
ss.x = s.x+2; v = v + 1i*LapDLP(s,ss,tau);
norm(u-v)
