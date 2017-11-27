function g = perispecinterp(f,N)
% PERISPECINTERP  resample periodically sampled function onto different grid.
%
% g = perispecinterp(f,N)
% inputs:  f - (row or column) vector length n (must be even) of input samples
%          N - desired output number of samples, must be even
% outputs: g - vector length N of interpolant (row or col, same size as f)
%
% Note on phasing: the output and input grid first entries align (ie, as if they
%  are both 0-indexed; note this matches setupquad).

% Barnett 6/27/16 renaming fftinterp from 9/5/14. Upsample case 11/16/17.
% todo: odd cases.

if nargin==0, test_perispecinterp; return; end
n = numel(f);
if N==n, g = f; return; end   % trivial case!
if mod(N,2)~=0 || mod(n,2)~=0, warning('Both N and n must be even; sorry'); end
F = fft(f(:).');    % row vector
if N>n              % upsample
  g = ifft([F(1:n/2) F(n/2+1)/2 zeros(1,N-n-1) F(n/2+1)/2 F(n/2+2:end)]);
else                % downsample
  g = ifft([F(1:N/2) F(end-N/2+1:end)]);
end
g = g*(N/n);   % factor from the ifft
if size(f,1)>size(f,2), g = g(:); end % make col vector
%%%%%%

function test_perispecinterp
n = 100;          % starting grid
x = 2*pi*(0:n-1)/n;
f = @(x) exp(sin(x));
Ns = [142 42];    % upsample test then downsample test
for i=1:2, N=Ns(i);
  g = perispecinterp(f(x),N);
  ge = f(2*pi*(0:N-1)/N);
  norm(g - ge)
end
