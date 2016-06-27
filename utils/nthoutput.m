function y = nthoutput(fun,n,varargin)
% NTHOUTPUT   return nth output of a function call
varargout = fun(varargin(:));
y = varargout(n);
