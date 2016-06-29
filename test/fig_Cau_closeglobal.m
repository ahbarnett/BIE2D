% convergence figs for Cauchy compensated evaluation for LSC2D paper Fig. 3.1.
% Barnett 10/28/13, repackaged 6/12/16. Path must include ../
%
% Also see self-test for Cau_closeglobal which prints a table of errors, no figs

a = .3; w = 5;         % smooth wobbly radial shape params
R = @(t) 1 + a*cos(w*t); Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);
s.inside = @(z) abs(z)<R(angle(z)); s.a = -0.1;  % interior pt far from bdry

side = 'e';   % either 'i' interior or 'e' exterior for which side to test

% holomorphic test func
%v = @(z) cos(z+1+0.5i); vp = @(z) -sin(z+1+0.5i);  % smoothish one
%v = @(z) cos(10*z+1+0.5i); vp = @(z) -10*sin(z+1+0.5i); % rapid osc + 1e5 size
%p=13; v = @(z) z.^p; vp = @(z) p*z.^(p-1); % monomial
%P=15; c = randn(1,P+1) + 1i*randn(1,P+1); v = @(z) polyval(c,z); % random poly
%vp = @(z) polyval((P:-1:1).*c(1:end-1),z);
%a = 2; v = @(z) exp(a*z); vp = @(z) a*exp(a*z);  % entire, interior
%a = 2; v = @(z) exp(a./z)-1; vp = @(z) -a./z.^2.*exp(a./z);  % entire in exterior

a = 1.1+1i; if side=='e', a = .1+.5i; end  % pole, dist 0.5 from G (.33 for ext)
v = @(z) 1./(z-a); vp = @(z) -1./(z-a).^2;   % used in paper.

% convergence at set of pts
z0 = s.Z(2*pi/4); % if N is multiple of 4, this is always a node
ds = logspace(0,-18,10)*(.1-1i); % displacements
if side=='e', ds = -ds; end % flip to outside
z = z0 + ds; z(end) = z0; % ray of pts heading to a node (make last hit exactly)
vz = v(z); vpz = vp(z); M = numel(z);

s = setupquad(s,200); %figure; plot(s.x,'k.-'); hold on; plot(z,'+'); axis equal
maxv = max(abs(v(s.x))); % size of v
fprintf('max |v| on bdry = %.3g, on test pts = %.3g\n',max(abs(v(s.x))), ...
        max(abs(vz)))
if 1, G = 1.3; if side=='e', G = 2.5; end % larger grid if outside
figure; g=-G:0.01:G; [xx yy]=meshgrid(g);zz=xx+1i*yy; title('Re v');
imagesc(g,g,real(v(zz))); hold on; plot(s.x,'k.-');
hold on; plot(z,'w+'); axis xy equal; caxis(maxv*[-1 1]); colorbar; end

o = []; %o.delta = 1e-2;  % no delta field uses default; 0 makes it non-bary
Ns = 20:20:200; err = nan(numel(Ns),numel(z)); errp = err; % N multiple of 4
for i=1:numel(Ns), N=Ns(i); s = setupquad(s,N);
  %d = repmat(s.x(:),[1 M])-repmat(z(:).',[N 1]); % displ mat
  %vc = sum(repmat(v(s.x).*s.cw,[1 M])./d,1)/(2i*pi); % naive Cauchy - so bad!
  [vc vcp] = Cau_closeglobal(z,s,v(s.x),side,o);       % my bary alg
  err(i,:) = vc.' - vz; errp(i,:) = vcp.' - vpz;
end

% FIGURE PLOTS

%s = quadr(s,200); figure; plot(s.x,'k.-'); hold on; plot(z,'+-'); axis equal;
%plot(a,'k*'); title('(a)'); axis off;
%set(gcf,'paperposition',[0 0 3 3.5]); print -depsc2 cau_int_geom.eps

maxv = 1; % not worth quibbling about

figure; imagesc(log10(abs(ds)),Ns,log10(abs(err/maxv))); colormap(goodbw);
caxis([-15 0]); axis xy;
xlabel('log_{10} dist(x,y_j)');ylabel('N');title('(b) log_{10} error v(x)');
%set(gcf,'paperposition',[0 0 2.7 3]); print -depsc2 cau_int_err.eps

figure; imagesc(log10(abs(ds)),Ns,log10(abs(errp/maxv))); colormap(goodbw);
caxis([-15 0]); colorbar;axis xy;
xlabel('log_{10} dist(x,y_j)');ylabel('N');title('(d) log_{10} error v''(x)');
%text(-4,200,'\delta=10^{-2}');
%set(gcf,'paperposition',[0 0 3.3 3]); print -depsc2 cau_int_errp.eps
