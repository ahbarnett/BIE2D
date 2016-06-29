function s = setupquad(s, N)
% SETUPQUAD  Set up periodic trapezoid quadrature & geom for smooth closed curve
%
% s = setupquad(s,N) where s is a struct containing a parametrization of the
%  curve in the form of function s.Z from [0,2pi) to the complex plane, uses
%  this to build the set of nodes, weights, speeds, curvatures, etc, using the
%  N-node PTR on [0,2pi).
%
% s = setupquad(s) where s contains at least the field s.x (column vector of
%  N node locations) generates all other fields by spectral differentiation.
%
% If s.Zp is present it is used as the function Z'; likewise s.Zpp is used for
%  Z''. The point of using these is that slightly more accurate normals and
%  curvatures may result compared to differentiation. On the other hand, Zp and
%  Zpp are not checked for plausibility, and are often complicated to code.
%
% Note: curves go counter-clockwise. Node j is at Z(2pi.(j-1)/N), ie 0-indexed.
%  FFT is used so is O(N log N), or O(N) if Z, Zp and Zpp available.
%
% Run without arguments, a self-test is done.
%
% Inputs:
%  s : struct containing either the field s.Z, a mapping from [0,2pi) to
%      the complex plane, which must be able to accept vectorized inputs,
%      or the field s.x containing column vector of the nodes.
%  N : number of nodes (not needed if s.x is present; however, if s.Z is also
%      present, N overrides the number of nodes in s.x).
% Outputs:
%  s : same struct with added column vector fields (nodes, weights, velocities,
%      curvatures, etc) needed for quadrature and Nystrom methods. Namely,
%      s.x nodes in complex plane, Z(s_j)
%      s.xp velocities Z'(s_j)
%      s.xp accelerations Z''(s_j)
%      s.t parameter values s_j for trapezoid rule, 2pi.(j-1)/N, j=1,...,N
%      s.nx outward unit normals
%      s.tang forward unit tangential vectors
%      s.sp speeds |Z'(s_j)|
%      s.w "speed weights" (2pi/N)*s.sp
%      s.cw velocity weights (ie complex speed)
%      s.cur curvatures kappa(s_j) (inverse bending radius)
%
% Example usage:
%
% s.Z = @(s) (1+0.3*cos(5*s).*exp(1i*s);                   % starfish param
% s = setupquad(s,100);
% figure; plot(s.x,'k.'); hold on; plot([s.x, s.x+0.2*s.nx].', 'b-'); axis equal
%
% Now check that normals from spectral differentiation are accurate:
%
% s.Zp = @(s) -1.5*sin(5*s).*exp(1i*s) + 1i*s.Z(s);        % Z' formula
% t = setupquad(s,100); norm(t.nx-s.nx)                 % should be small
%
% Also see: PERISPECDIFF.

% (c) Alex Barnett 10/8/14, name changed to avoid conflict w/ mpspack 6/12/16.
% 0-indexed to match interp, 6/29/16

if nargin==0, test_setupquad; return; end
if nargin>1           % use N from args
  s.t = (0:N-1)'*(2*pi/N);
  if isfield(s,'Z'), s.x = s.Z(s.t); end    % use formula
  if N~=length(s.x), error('N differs from length of s.x; that sucks!'); end 
elseif isfield(s,'x')
  s.x = s.x(:);          % ensure col vec
  N = length(s.x);
  s.t = (0:N-1)'*(2*pi/N); % we don't know the actual params, but choose this
else
  error('Need to provide at least s.Z and N, or s.x. Neither found!');
end
if isfield(s,'Zp'), s.xp = s.Zp(s.t); else, s.xp = perispecdiff(s.x); end
if isfield(s,'Zpp'), s.xpp = s.Zpp(s.t); else, s.xpp = perispecdiff(s.xp); end
% Now local stuff that derives from x, xp, xpp at each node...
s.sp = abs(s.xp);
s.tang = s.xp./s.sp;
s.nx = -1i*s.tang;
s.cur = -real(conj(s.xpp).*s.nx)./s.sp.^2;  % recall real(conj(a)*b) = "a dot b"
s.w = (2*pi/N)*s.sp;
s.cw = (2*pi/N)*s.xp;  % complex weights (incl complex speed)

%%%%%%%%%%%%%%%%%%%%%%%%%
function test_setupquad                % not very extensive! Barnett 10/8/14
Z = @(s) (1+0.3*cos(5*s)).*exp(1i*s);                   % starfish param
s.Z = Z;
n = 100;
s = setupquad(s,n);
s = []; s.x = Z((0:n-1)/n*2*pi); s = setupquad(s);  % testing s.x input only
figure; plot(s.x,'k.'); hold on; plot([s.x, s.x+0.2*s.nx].', 'b-'); axis equal
% Now check that normals from spectral differentiation are accurate:
Zp = @(s) -1.5*sin(5*s).*exp(1i*s) + 1i*Z(s);        % Z' formula
s.Zp = Zp;
t = setupquad(s); norm(t.nx-s.nx)             % should be small
s = []; s.x = 3; s.Z = Z; s = setupquad(s,100);    % N should override s.x
%s = []; s.x = Z((1:50)'/50*2*pi); s=setupquad(s,100); should fail
