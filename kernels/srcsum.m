function [A B] = srcsum(kernel, trlist, t, s, varargin)
% SRCSUM   Sum any kernel evaluation matrix over a set of source translations
%
% [A B] = srcsum(kernel, trlist, t, s)
% [A B] = srcsum(kernel, trlist, t, s, param)
%
% kernel: function of form [A ...] = kernel(t,s) or kernel(t,s,param).
% trlist = list of complex numbers giving source translations
% t = target segment struct
% s = source segment struct

% Barnett 6/12/16

xo = s.x; s.x = xo + trlist(1);
if nargout==1
  A = kernel(t,s,varargin{:});  % since we don't know the size can't prealloc
  for i=2:numel(trlist)
    s.x = xo + trlist(i);
    A = A + kernel(t,s,varargin{:});            % todo: phases for Helm
  end
else
  [A B] = kernel(t,s,varargin{:});
  for i=2:numel(trlist)
    s.x = xo + trlist(i);
    [Ai Bi] = kernel(t,s,varargin{:});
    A = A + Ai; B = B + Bi;
  end
end
