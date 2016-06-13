function g = perispecdiff(f)
% PERISPECDIFF - use FFT to take periodic spectral differentiation of vector
%
% g = perispecdiff(f) returns g the derivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points). Can be row or col vec, and output is same shape.
%
% Without arguments, does a self-test.

% Barnett 2/18/14
if nargin==0, test_perispecdiff; return; end
N = numel(f);
if mod(N,2)==0   % even
  g = ifft(fft(f(:)).*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].');
else
  g = ifft(fft(f(:)).*[0 1i*(1:(N-1)/2) 1i*((1-N)/2:-1)].');
end
g = reshape(g,size(f));

%%%%%%
function test_perispecdiff
N = 50; tj = 2*pi/N*(1:N)';
f = sin(3*tj); fp = 3*cos(3*tj);   % trial periodic function & its deriv
norm(fp-perispecdiff(f))
