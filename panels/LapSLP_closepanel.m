function [A, A1, A2] = LapSLP_closepanel(t,s,a,b,side)
% LAPSLP_CLOSEPANEL - SLP val+grad close-eval Helsing "special quadrature" matrix
%
% [A] = LapSLP_closepanel(t,s,a,b) returns
%  returns numel(t.x)-by-numel(s.x) matrix which maps SLP values at the nodes
%  s.x to potential at the targets t.x, both given as lists of
%  points in the complex plane.
%  The matrix is the quadrature approximation to evaluation of
%    (1/(2*pi)) * int_Gamma log |x-y| sigma(y) dy,
%  ie the 2D Laplace SLP.
%
% [A An] = LapSLP_closepanel(t,s,a,b) also gives target normal-derivs (needs t.nx)
%  *** TODO. UNTESTED DERIVS  ***
% [A A1 A2] = LapSLP_closepanel(t,s,a,b) also gives target x,y-derivs
%
% Inputs: t = target seg struct (with column-vec t.x targets in complex plane)
%         s = src node seg struct (with s.x, s.wxp, s.nx)
%         a = panel start, b = panel end, in complex plane.
%         side = 'i' (interior) or 'e' (exterior). [optional]
% Output: A (n_targ * n_src) is source-to-target value matrix
%         An or A1, A2 = source to target normal-deriv (or x,y-deriv) matrices
%
% For test: see ../test/testGRFLap_panels.m for now

% Notes:
% 0) adapted from Dspecialquad.m in stokes-panel-quad, orig laplacecloseeval.m
% 1) Efficient only if multiple targs, since O(p^3).
% 2) See Helsing-Ojala 2008 (special quadr Sec 5.1-2),
%  Helsing 2009 mixed (p=16), and Helsing's tutorial demo11b.m M1IcompRecFS().
% 3) real part is taken, which prevents the Stokes extension using complex tau.
% 4) not tidy. still uses curved decision boundary (gam), to clear up.
%
% Authors: Alex Barnett (2013-2021), based on Johan Helsing.
% tweaks by Bowei Wu, Hai Zhu.

if nargin<5, side = 'i'; end     % interior or exterior
zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc;  % transformed src nodes, targ pts
N = numel(x);                            % # of targets
p = numel(s.x);                          % assume panel order is # nodes
c = (1-(-1).^(1:p))./(1:p);              % Helsing c_k, k = 1..p.
V = ones(p,p); for k=2:p, V(:,k) = V(:,k-1).*y; end  % Vandermonde mat @ nodes
P = zeros(p+1,N);      % Build P, Helsing's p_k vectorized on all targs...
d = 1.1; inr = abs(x)<=d; ifr = abs(x)>d;      % near & far treat separately
% compute P up to p+1 instead of p as in DLP, since q_k needs them:
%gam = 1i;  % original choice: branch cut is semicircle behind panel
gam = exp(1i*pi/4);  % smaller makes cut closer to panel. barnett 4/17/18
% note gam = 1 fails, and gam = -1 put cuts in the domain.
if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
P(1,:) = log(gam) + log((1-x)./(gam*(-1-x)));  % init p_1 for all targs int

% upwards recurrence for near targets, faster + more acc than quadr...
% (note rotation of cut in log to -Im; so cut in x space is lower unit circle)
Nn =  numel(find(inr));
if Nn ~= 0  % Criterion added by Hai Zhu 08/24/16 to ensure inr not empty
    for k=1:p
        P(k+1,inr) = x(inr).'.*P(k,inr) + c(k); 
    end  % recursion for p_k
end
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf =  numel(find(ifr)); wxp = s.wxp/zsc; % rescaled complex speed weights
if Nf>0 % Criterion added by Bowei Wu 03/05/15 to ensure ifr not empty
    P(end,ifr) = sum(((wxp.*(V(:,end).*y(:)))*ones(1,Nf))./bsxfun(@minus,y,x(ifr).'));  % int y^p/(y-x)
    for ii = p:-1:2
        P( ii,ifr) = (P(ii+1,ifr)-c(ii))./x(ifr).'; % backward recursion
    end
end

Q = zeros(p,N); % compute q_k from p_k via Helsing 2009 eqn (18)... (p even!)
% Note a rot ang appears here too...  4/17/18
%gam = exp(1i*pi/4); % 1i;  % moves a branch arc as in p_1
%if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
Q(1:2:end,:) = P(2:2:end,:) - repmat(log((1-x.').*(-1-x.')),[ceil(p/2) 1]); % guessed!
% (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
Q(2:2:end,:) = P(3:2:end,:) - repmat(log(gam) + log((1-x.')./(gam*(-1-x.'))),[floor(p/2) 1]);  % same cut as for p_1
% Seems like abs fails - we must be using complex SLP ? :
%Q(1:2:end,:) = P(2:2:end,:) - repmat(log(abs(1-x.'))+log(abs(-1-x.')),[p/2 1]);
% (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
%Q(2:2:end,:) = P(3:2:end,:) - repmat(log(abs(1-x.'))-log(abs(-1-x.')),[p/2 1]);
Q = Q.*repmat(1./(1:p)',[1 N]); % k even
warning('off','MATLAB:nearlySingularMatrix'); % solve for special weights...
A = real((V.'\Q).'.*repmat((1i*s.nx)',[N 1])*zsc)/(2*pi*abs(zsc));
A = A*abs(zsc) - log(abs(zsc))/(2*pi)*repmat(abs(s.wxp)',[N 1]); % unscale, yuk
if nargout>1
    P = P(1:end-1,:);  % trim P back to p rows since kernel is like DLP
    Az = (V.'\P).'*(1/(2*pi)).*repmat((1i*s.nx)',[N 1]); % solve spec wei
    if nargout == 2
        A1 = Az;
    else
        A1 = real(Az); A2 = -imag(Az);
    end
end
end

