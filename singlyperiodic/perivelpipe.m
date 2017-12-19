function perivelpipe
% singly-periodic in 2D velocity-BC Stokes BVP, ie "pipe" geom w/ press drop
% Rewrite of 2014 codes using BIE2D, as in DPLS codes. Dense E fill for now,
% DLP only on U,D (since wrapping S quad would be harder).  12/17/17

% Issues:
% * we don't have close eval for open (periodized?) segments yet.
% * rewrite using periodized Stokes FMM, which should have a direct option too.

clear; v=1;  % verbosity=0,1,2,3
expt='t';    % expt='t' test known soln, 'd' driven no-slip demo

% set up upper and lower walls
uc.e1 = 2*pi;   % unitcell, e1=lattice vector as a complex number
N = 70;     % pts per top and bottom wall (enough for 1e-15 in this case)
U.Z = @(t) t + 1i*(1.5+sin(t)); U.Zp = @(t) 1 + 1i*cos(t);
U.Zpp = @(t) -1i*sin(t);
U = setupquad(U,N);   % t=0 is left end, t=2pi right end
U.nx = -U.nx; U.cur = -U.cur; U.xp = -U.xp; U.tang=-U.tang; U.cw=-U.cw; % correct for sense of U, opp from periodicdirpipe
% reorder pts, makes no diff...
%U.x=U.x(end:-1:1); U.nx=U.nx(end:-1:1); U.xp=U.xp(end:-1:1); U.xpp=U.xpp(end:-1:1); U.tang=U.tang(end:-1:1); U.cw=U.cw(end:-1:1); U.w=U.w(end:-1:1); U.sp=U.sp(end:-1:1); U.cur=U.cur(end:-1:1); 

D.Z = @(t) t + 1i*(-1.5+cos(2*t)); D.Zp = @(t) 1 - 2i*sin(2*t);
D.Zpp = @(t) - 4i*cos(2*t);
D = setupquad(D,N);   % same dir (left to right); pipe wall normals outwards
inside = @(z) imag(z-D.Z(real(z)))>0 & imag(z-U.Z(real(z)))<0; % ok U,D graphs
s = mergesegquads([U,D]);
zt = 2+1i;    % point to test u soln at (expt='t' only)

% set up left and right walls
uc.nei = 1; % how many nei copies either side (use eg 1e3 to test A w/o AP)
uc.trlist = uc.e1*(-uc.nei:uc.nei);  % list of translations for direct images
m = 20;    % pts per side wall
[x w] = gauss(m); x = (1+x)/2; w = w'/2; % quadr on [0,1]
H = U.Z(0)-D.Z(0); L.x = D.Z(0) + H*x; L.nx = 0*L.x+1; L.w = H*w; % left side
R = L; R.x = L.x+uc.e1; % right side
uc.L = L; uc.R = R;

% set up aux periodizing basis
proxyrep = @StoSLP;      % sets proxy pt type via a kernel function call
Rp = 1.1*2*pi; M = 70;    % # proxy pts (2 force comps per pt, so 2M dofs)
p.x = pi + Rp*exp(2i*pi*(0:M-1)'/M); p = setupquad(p);     % proxy pts
if v>1, figure; showsegment({U,D},uc.trlist); showsegment({L,R}); plot(p.x,'r+'); plot(zt,'go'); end

mu = 0.7;                                           % fluid viscosity
if expt=='t' % Exact soln: either periodic or plus fixed pressure drop / period:
  %ue = @(x) [1+0*x; -2+0*x]; pe = @(x) 0*x; % the exact soln: uniform rigid flow, constant pressure everywhere (no drop)
  h=.2; ue = @(x) h*[imag(x).^2;0*x]; pe = @(x) h*2*mu*real(x); % horiz Poisseuil flow (has pres drop)
  disp('expt=t: running known Poisseuil flow BVP...')
  vrhs = ue(s.x);         % bdry vel data: NB ordering Ux,Dx,Uy,Dy !
  jump = pe(uc.e1)-pe(0); % known pressure growth across one period (a number)
elseif expt=='d'
  vrhs = zeros(4*N,1);     % no-slip BCs, ie homog vel data on U,D
  jump = -1;              % given pressure driving, for flow +x (a number)
  disp('expt=d: solving no-slip pressure-driven flow in pipe...')
end

tic
E = ELSmatrix(s,p,proxyrep,mu,uc);                % fill
if v>2, S = svd(E); disp('last few sing vals of E:'), S(end-5:end)
  figure; imagesc(E); title('E'); colorbar; end
%sum(rhs(1:N).*real(U.nx)), sum(rhs(N+1:2*N).*imag(U.nx)) % fluid cons thru U
%sum(rhs(2*N+1:3*N).*real(D.nx)), sum(rhs(3*N+1:end).*imag(D.nx)) % fluid cons thru D
Tjump = -jump * [real(R.nx);imag(R.nx)]; % traction driving growth (vector func)
erhs = [vrhs; zeros(2*m,1);Tjump];       % expanded lin sys RHS
warning('off','MATLAB:nearlySingularMatrix')  % backward-stable ill-cond is ok!
warning('off','MATLAB:rankDeficientMatrix')
lso.RECT = true;  % linsolve opts, forces QR even when square
co = linsolve(E,erhs,lso);                           % direct bkw stable solve
toc
fprintf('resid norm = %.3g\n',norm(E*co - erhs))
sig = co(1:4*N); psi = co(4*N+1:end);
fprintf('density norm = %.3g, proxy norm = %.3g\n',norm(sig), norm(psi))
[ut pt] = evalsol(s,p,proxyrep,mu,uc,zt,co);
if expt=='t'                         % check vs known soln at the point zt
  fprintf('u velocity err at zt = %.3g\n', norm(ut-ue(zt)))
else, fprintf('u velocity at zt = [%.15g, %.15g]\n', ut(1),ut(2))
end

if v   % plots
  nx = 60; gx = 2*pi*((1:nx)-0.5)/nx; ny = nx; gy = gx - pi; % plotting grid
  [xx yy] = meshgrid(gx,gy); t.x = xx(:)+1i*yy(:); Mt = numel(t.x);
  di = reshape(inside(t.x),size(xx));  % boolean if inside domain
  if expt=='t', ueg = ue(t.x);      % evaluate & show known soln on grid...
    ue1 = reshape(ueg(1:Mt),size(xx)); ue2 = reshape(ueg(Mt+(1:Mt)),size(xx));
    peg = reshape(pe(t.x),size(xx));
    figure; imagesc(gx,gy, peg); colormap(jet(256)); colorbar;
    hold on; quiver(gx,gy, ue1,ue2, 10);
  else, figure; end
  showsegment({U,D},uc.trlist); showsegment({L,R}); plot(p.x,'r+'); plot(zt,'go');
  text(-.5,0,'L');text(2*pi+.5,0,'R');text(4,-1,'D');text(pi,0.5,'U');
  title('geom'); if expt=='t',title('geom and (u,p) known soln'); end
  % eval and plot soln...
  ug = nan(size([t.x;t.x])); pg = nan(size(t.x)); ii = inside(t.x);
  [ug([ii;ii]) pg(ii)] = evalsol(s,p,proxyrep,mu,uc,t.x(ii),co);
  u1 = reshape(ug(1:Mt),size(xx)); u2 = reshape(ug(Mt+(1:Mt)),size(xx));
  pp = reshape(pg,size(xx)); pp = pp - pp(ceil(ny/2),1); % zero p mid left edge
  figure; imagesc(gx,gy, pp); colormap(jet(256));
  caxis(sort([0 jump])); colorbar;
  hold on; quiver(gx,gy, u1,u2, 10); title('soln (u,p), w/o close-eval scheme')
  showsegment({U,D}); showsegment({L,R}); plot(zt,'go');
  if expt=='t'     % show error vs known...
    i=ceil(ny/2); j=ceil(nx/4); % index interior pt to get pres const
    pp = pp - pp(i,j) + peg(i,j);     % shift const in pres to match known
    eg2 = sum([u1(:)-ue1(:),u2(:)-ue2(:)].^2,2); % squared ptwise vector L2 errs
    figure; subplot(1,2,1); imagesc(gx,gy,log10(reshape(eg2,size(xx)))/2);
    axis xy equal tight; caxis([-16 0]); colorbar; hold on; plot(U.x,'k.-');
    plot(D.x,'k.-'); title('peri vel Stokes BVP: log_{10} u err')
    subplot(1,2,2); imagesc(gx,gy,log10(abs(pp-reshape(peg,size(xx)))));
    axis xy equal tight; caxis([-16 0]); colorbar; hold on;
    plot(U.x,'k.-'); plot(D.x,'k.-'); title('log_{10} p err (up to const)');
  end
end

%keyboard % don't forget to use dbquit to finish otherwise trapped in debug mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u p] = evalsol(s,pr,proxyrep,mu,U,z,co,close) % eval soln rep u,p
% z = list of targets as C values. pr = proxy struct, s = source curve struct
% co = full coeff vec, U = unit cell struct, mu = viscosity. DLP only.
% Note: not for self-eval since srcsum2 used, 6/30/16. Taken from fig_stoconvK1.
% TODO: implement close eval.
if nargin<8, close=0; end              % default is plain native quadr
if close, error('close not implemented for open segments!'); end
z = struct('x',z);                     % make targets z a segment struct
N = numel(s.x);
sig = co(1:2*N); psi = co(2*N+1:end);  % split into sig (density) & psi (proxy)
if close
  D = @(t,s,mu,dens) StoDLP_closeglobal(t,s,mu,dens,'e'); % exterior
else, D = @StoDLP; end                 % NB now native & close have same 4 args
if nargout==1                          % don't want pressure output
  u = proxyrep(z,pr,mu,psi);           % init sol w/ proxies (always far)
  u = u + srcsum2(D,U.trlist,[],z,s,mu,sig);
else  
  [u p] = proxyrep(z,pr,mu,psi);       % init sol w/ proxies (always far)
  [uD pD] = srcsum2(D,U.trlist,[],z,s,mu,sig);
  u = u + uD;
  p = p + pD;
end

%function J = evalflux(s,p,proxyrep,mu,sd,uc,co)      % TODO
% direct L wall integral for now, needs close eval, which we don't have...
% 12/17/17
%close=0; u = evalsol(s,p,proxyrep,mu,sd,uc,uc.L.x,co,close);
%J = sum(u.*uc.L.w(:));

function [E A B C Q] = ELSmatrix(s,p,proxyrep,mu,uc)
% builds matrix blocks for Stokes extended linear system, D rep w/ Kress self
N = numel(s.x);
A = -eye(2*N)/2 + srcsum(@StoDLP,uc.trlist,[],s,s,mu);  % notes: DLP gives int JR term; srcsum self is ok
B = proxyrep(s,p,mu);     % map from proxy density to vel on curve
d = uc.e1*uc.nei;         % src transl to use
[CLD,~,TLD] = srcsum(@StoDLP,d,[],uc.L,s,mu);
[CRD,~,TRD] = srcsum(@StoDLP,-d,[],uc.R,s,mu);
C = [CRD-CLD; TRD-TLD];
[QL,~,QLt] = proxyrep(uc.L,p,mu); [QR,~,QRt] = proxyrep(uc.R,p,mu); % vel, tract
Q = [QR-QL; QRt-QLt];
E = [A B; C Q];
