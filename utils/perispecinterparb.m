function [g gfun] = perispecinterparb(f,t)
% PERISPECINTERPARB  interpolate grid-sampled periodic function to arb targets
%
% [g gfun] = perispecinterparb(f,t)
% inputs:  f - (row or column) vector length n (must be even) of input samples
%          t - list of ordinates to evaluate at (periodic wrt [0,2pi))
% outputs: g - values of interpolant at each t (row or col, whichever t was)
%          gfun - function handle which can map future t values to interpolant
%                 In particular, g = gfun(t) holds.
%
% Example: see self-test, which is done when called without arguments
%
% Note on phasing: the input grid first entry is at t=0, as if 0-indexed
% (this matches setupquad).
% Note on speed: naive complex summation used, so cost is O(numel(t).numel(f))

% Todo: odd n case
% Barnett 2/25/21
if nargin==0, test_perispecinterparb; return; end

n = numel(f);
if mod(n,2)~=0, warning('n must be even; sorry'); end
fhat = fft(f(:).')/n;    % row vector of Fourier series coeffs
gfun = @(t) fourierseries(fhat,t);
g = gfun(t);

function z = fourierseries(zhat,t)   % eval cmplx F series (t any shape matrix)
% *** todo: block out if t*n is too large, to keep in fast L2 cache?
sizet = size(t); t = t(:); zhat = zhat(:);                     % col vecs
n = numel(zhat);                 % must be even for now
z = 0*t + zhat(1);               % allocate and zero mode
% use outer products in the trig args... and zero-pad for k=-n/2
z = z + cos(t*(1:n/2)) * (zhat(2:n/2+1) + [zhat(end:-1:n/2+2);0]);
z = z + sin(t*(1:n/2)) * (1i*zhat(2:n/2+1) - 1i*[zhat(end:-1:n/2+2);0]);
z = reshape(z,sizet);

%%%%%%
function test_perispecinterparb
n = 70;          % input grid: for below f, expect err ~ exp(-striphalfwid*n/2)
                 % See R. Kress, Linear Integral Equations (1989), Thm. 11.5.
x = 2*pi*(0:n-1)/n;
f = @(x) cot((x+2.0+1i)/2); % generic complex 2pi-per func, holom striphalfwid=1
t = 2*pi*rand(20,10);       % targets
[g gfun] = perispecinterparb(f(x),t);
disp(max(abs(g - f(t))) / max(abs(g)))   % rel max err, cf exp(-1.0*n/2)

% timing
n=300; [~,gfun] = perispecinterparb(f(2*pi*(0:n-1)/n),t);
m = 1e4; r = 10; t=rand(m,1); tic; for i=1:r; g=gfun(t); end, t=toc;
fprintf('%.3g mode-evals/sec\n',n*m*r/t)    % get 1e8, big or small cases
