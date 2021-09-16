function [varargout] = LapDLP_closepanel(t,s,a,b,side, meth)
% LAPDLP_CLOSEPANEL - real DLP close-eval matrix (Helsing or sing-swap)
%
% [A] = LapDLP_closepanel(t,s,a,b)
%  returns numel(t.x)-by-numel(s.x) matrix which maps DLP values at the nodes
%  s.x to potential at the targets t.x, both given as lists of
%  points in the complex plane.
%  The matrix is the quadrature approximation to evaluation of
%    Re (1/(2i*pi)) * int_Gamma sigma(y)/(y-x) dy,
%  ie the 2D Laplace DLP acting on real densities (not the Cauchy complex val).
%
% [A] = LapDLP_closepanel(t,s,a,b,side) allows exterior if side='e', otherwise
%  interior side='i' is assumed.
%
% [A] = LapDLP_closepanel(t,s,a,b,side,meth) also controls method:
%   'h' = Helsing-Ojala (matches complex polynomial on curve)
%   's' = singulary-swap of af Klinteberg-Barnett (matches on flat panel)
%
% [A Az] = LapDLP_closepanel(t,s,a,b) also gives target complex gradient
%   *** untested since converted to Re part *** TO FIX
% [A Az Azz] = LapDLP_closepanel(t,s,a,b) also gives target gradient and Hessian (needs t.nx)
% [A A1 A2 A3 A4] = LapDLP_closepanel(t,s,a,b) also gives target x,y-derivs
%                   grad_t(D) = [A1, A2]; hess_t(D) = [A3, A4; A4, -A3];
%
% Inputs: t = target seg struct (with column-vec t.x targets in complex plane)
%         s = src node seg struct (with s.x, s.w, s.wxp; amazingly, s.nx not used!)
%         a = panel start, b = panel end, in complex plane.
%         [optional] side = 'i' (interior) or 'e' (exterior)
%         [optional] meth = 'h' (Helsing-Ojala, default)
%                           's' (af Klinteberg-Barnett, singularity-swap)
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
% 4) not tidy. still uses arc decision boundary (gam). to clear up.
%
% Authors: Alex Barnett (2013-2021); tweaks to HO code by Bowei Wu, Hai Zhu.

if nargin<5, side = 'i'; end     % interior or exterior
if nargin<6, meth = 'h'; end
% wrap the below routines...
if meth=='h'
  varargout{1:nargout} = LapDLP_closepanel_HO(t,s,a,b,side);
elseif meth=='s'
  varargout{1:nargout} = LapDLP_closepanel_SS(t,s,a,b,side);
end
%%%%%%%%%
  
function varargout = LapDLP_closepanel_HO(t,s,a,b,side)    % Re wrapper
% true DLP, taking real part.
varargout{1:nargout} = LapDLPcmplx_closepanel_HO(t,s,a,b,side);
for i=1:nargout,
  varargout{i} = real(varargout{i});
end

function [A, A1, A2, A3, A4] = LapDLPcmplx_closepanel_HO(t,s,a,b,side)
% the Helsing-Ojala version; see docs above. Complex (no taking Re part).
zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc;  % transformed src nodes, targ pts
%figure; plot(x,'.'); hold on; plot(y,'+-'); plot([-1 1],[0 0],'ro'); % debug
N = numel(x);                            % # of targets
p = numel(s.x);                          % assume panel order is # nodes
if N*p==0
    A = 0; A1=0; A2=0;
    return
end
c = (1-(-1).^(1:p))./(1:p);              % Helsing c_k, k = 1..p.
V = ones(p,p); for k=2:p, V(:,k) = V(:,k-1).*y; end  % Vandermonde mat @ nodes
P = zeros(p,N);      % Build P, Helsing's p_k vectorized on all targs...
d = 1.1; inr = abs(x)<=d; ifr = abs(x)>d;      % near & far treat separately
%gam = 1i;
gam = exp(1i*pi/4);  % smaller makes cut closer to panel. barnett 4/17/18
if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
P(1,:) = log(gam) + log((1-x)./(gam*(-1-x)));  % init p_1 for all targs int

% upwards recurrence for near targets, faster + more acc than quadr...
% (note rotation of cut in log to -Im; so cut in x space is lower unit circle)
if N>1 || (N==1 && inr==1) % Criterion added by Bowei Wu 03/05/15 to ensure inr not empty
    for k=1:p-1, P(k+1,inr) = x(inr).'.*P(k,inr) + c(k); end  % recursion for p_k
end
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf = numel(find(ifr)); wxp = s.wxp/zsc;   % rescaled complex speed weights

if Nf>0 % Criterion added by Bowei Wu 03/05/15 to ensure ifr not empty
    P(end,ifr) = sum(((wxp.*V(:,end))*ones(1,Nf))./bsxfun(@minus,y,x(ifr).'));
    for ii = p-1:-1:2
        P( ii,ifr) = (P(ii+1,ifr)-c(ii))./x(ifr).';
    end
end

warning('off','MATLAB:nearlySingularMatrix');
%A = real((V.'\P).'*(1i/(2*pi)));         % solve for special quadr weights
A = ((V.'\P).'*(1i/(2*pi)));         % do not take real for the eval of Stokes DLP non-laplace term. Bowei 10/19/14
%A = (P.'*inv(V))*(1i/(2*pi));   % equiv in exact arith, but not bkw stable.
if nargout>1
    R =  -(kron(ones(p,1),1./(1-x.')) + kron((-1).^(0:p-1).',1./(1+x.'))) +...
        repmat((0:p-1)',[1 N]).*[zeros(1,N); P(1:p-1,:)];  % hypersingular kernel weights of Helsing 2009 eqn (14)
    Az = (V.'\R).'*(1i/(2*pi*zsc));  % solve for targ complex-deriv mat & rescale
    A1 = Az;
    if nargout > 2
        S = -(kron(ones(p,1),1./(1-x.').^2) - kron((-1).^(0:p-1).',1./(1+x.').^2))/2 +...
            repmat((0:p-1)',[1 N]).*[zeros(1,N); R(1:p-1,:)]/2; % supersingular kernel weights
        Azz = (V.'\S).'*(1i/(2*pi*zsc^2));
        if nargout > 3
            A1 = real(Az); A2 = -imag(Az);  % note sign for y-deriv from C-deriv
            A3 = real(Azz); A4 = -imag(Azz);    
        else
            A1 = Az; A2 = Azz; 
        end
    end
end


function [A, A1, A2, A3, A4] = LapDLP_closepanel_SS(t,s,a,b,side)
% the singularity-swap version; see docs above.
% Just the nodes (and stdpan nodes) are needed, no analytic parameterization.

% todo: vectorize over targs t
p=numel(s.x);
[stdpan.t, stdpan.w] = gauss(p);   % todo: pass this in or look up
stdpan.x = stdpan.t; stdpan.wxp = stdpan.w';    % bare min for stdpan struct
[t0.x, c, cp] = panel_preimage(s.x,stdpan.x,t.x);  % t0=param-plane targ struct
A = LapDLPcmplx_closepanel_HO(t0,stdpan,-1,1,side);  % Helsing eval mat flat panel

ypj = polyval(cp, stdpan.t);  % get map derivs at nodes
% right-diag-scale the (complex) HO matrix to get the SS matrix...
scvec = (stdpan.t-t0.x)./(s.x-t.x).*ypj;
A = A*diag(scvec);
A = real(A);

% *** todo: deriv matrices
