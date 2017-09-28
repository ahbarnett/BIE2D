function g = perispecint(f)
% PERISPECINT - use FFT to take periodic spectral antiderivative of vector
%
% g = perispecint(f) returns g an antiderivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points). Can be row or col vec, and output is same shape.
%  If sum(f)=0 then g is smoothly periodic; otherwise, there is a sawtooth
%  jump in g. The offset of g is arbitrary.
%
% Without arguments, does a self-test.
%
% Also see: PERISPECDIFF

% Barnett 9/28/17
if nargin==0, test_perispecint; return; end
N = numel(f);
fbar = mean(f); f = f-fbar;
if mod(N,2)==0   % even
  g = ifft(fft(f(:)).*[0 1./(1i*(1:N/2-1)) 0 1./(1i*(-N/2+1:-1))].');
else
  g = ifft(fft(f(:)).*[0 1./(1i*(1:(N-1)/2)) 1./(1i*((1-N)/2:-1))].');
end
g = g + (1:N)'*(fbar*2*pi/N);  % add a sawtooth (could add arbitrary offset too)
g = reshape(g,size(f));

%%%%%%
function test_perispecint
N = 50; tj = 2*pi/N*(1:N)';
f = sin(3*tj); fp = 3*cos(3*tj);   % trial periodic function & its deriv
% since offset arbitrary, must measure and subtract...
F = perispecint(fp); off = F(1)-f(1); norm(F-off-f)   % zero-mean case for fp
f = sin(3*tj)+tj; fp = 3*cos(3*tj)+1;   % trial periodic function & its deriv
F = perispecint(fp); off = F(1)-f(1); norm(F-off-f)   % nonzero-mean case for fp
