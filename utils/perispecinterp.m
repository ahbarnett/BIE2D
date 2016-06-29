function g = perispecinterp(f,N)
% PERISPECINTERP   resample periodically sampled function on finer uniform grid.
%
% g = perispecinterp(f,N)
% inputs:  f - (row or column) vector length n (must be even) of samples
%          N - desired output number of samples, must be >= n and even
% outputs: g - vector length N of interpolant, (row or col as f was)
%
% Note on phasing: the output and input grid first entry align (ie, as if they
%  are both 0-indexed; note this matches setupquad)

% Barnett 6/27/16 renaming fftinterp from 9/5/14
% todo: * downsample case N<n.  * odd cases.

if nargin==0, test_perispecinterp; return; end
n = numel(f);
if N==n, g = f; return; end
if mod(N,2)~=0 || mod(n,2)~=0, warning('N and n must be even, sorry'); end
F = fft(f(:).');    % row vector
g = ifft([F(1:n/2) F(n/2+1)/2 zeros(1,N-n-1) F(n/2+1)/2 F(n/2+2:end)]);
g = g*(N/n);   % factor from the ifft
if size(f,1)>size(f,2), g = g(:); end % make col vector
%%%%%%

function test_perispecinterp
n = 50;
N = 100;
x = 2*pi*(0:n-1)/n;
f = @(x) exp(sin(x));
g = perispecinterp(f(x),N);
ge = f(2*pi*(0:N-1)/N);
% g ./ ge
norm(g - ge)
