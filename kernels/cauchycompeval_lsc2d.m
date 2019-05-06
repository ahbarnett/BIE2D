function [vc, vcp] = cauchycompeval_lsc2d(x,s,vb,side,o)
% cauchycompeval_lsc2d.   Legacy variant of cauchycompeval, from [lsc2d] paper.
%
% See cauchcompeval
%
% The only new input needed is s.a, an interior point far from the boundary.
%
% The algorithm is that of Ioakimidis et al BIT 1991 for interior, and,
% for the exterior, a modified version using 1/(z-a) in place of the function 1.
% For the derivative, the formula is mathematically the derivative of the
% barycentric formula of Schneider-Werner 1986 described in Berrut et al 2005,
% but using the complex quadrature weights instead of the barycentric weights.
% However, since this derivative is not itself a true barycentric form, a hack
% is needed to compute the difference v_j-v(x) in a form where roundoff error
% cancels correctly. The exterior case is slightly more intricate. When
% a target coincides with a node j, the true value v_j (or Schneider-Werner
% formula for v'(x_j)) is used.
%
% Alex Barnett 10/22/13 based on cauchycompevalint, pole code in lapDevalclose.m
% 10/23/13 node-targ coincidences fixed, exterior true barycentric discovered.

% todo: check v' ext at the node (S-W form), seems to be wrong.

if nargin<5, o = []; end
if size(vb,2)>1              % matrix versions
  if nargout==1, vc = cauchycompmat_lsc2d(x,s,vb,side,o);
  else, [vc vcp] = cauchycompmat_lsc2d(x,s,vb,side,o); end
  return
end
if ~isfield(o,'delta'), o.delta = 1e-2; end  % currently won't ever be run.
% delta = dist param where $O(N^2.M) deriv bary form switched on.
% Roughly deriv errors are then limited to emach/mindist
% The only reason to decrease mindist is if lots of v close
% nodes on the curve with lots of close target points.

cw = s.cw;                    % complex speed weights
N = numel(s.x); M = numel(x);

if nargout==1  % no deriv wanted... (note sum along 1-axis faster than 2-axis)
  comp = repmat(cw(:), [1 M]) ./ (repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1]));
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  I0 = sum(repmat(vb(:),[1 M]).*comp); J0 = sum(pcomp); % Ioakimidis notation
  vc = I0./J0;                        % bary form
  if side=='e', vc = vc./(x(:).'-s.a); end         % correct w/ pole
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), vc(ii(l)) = vb(jj(l)); end % replace each hit w/ corresp vb
  
else           % 1st deriv wanted...
  invd = 1./(repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1])); % 1/displacement mat
  comp = repmat(cw(:), [1 M]) .* invd;
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  I0 = sum(repmat(vb(:),[1 M]).*comp); J0 = sum(pcomp);
  if side=='e', prefac = 1./(x(:).'- s.a); else prefac = 1; end
  vc = prefac .* I0./J0;                        % bary form (poss overall pole)
  dv = repmat(vb(:),[1 M]) - repmat(vc(:).',[N 1]); % v value diff mat
  [jj ii] = ind2sub(size(invd),find(abs(invd) > 1/o.delta)); % bad pairs indices
  if side=='e' % exterior:
    for l=1:numel(jj), j=jj(l); i=ii(l); % loop over node-targ pairs too close
      p = sum(comp(:,i).*(vb(j)./(s.x(:)-s.a)-vb(:)./(s.x(j)-s.a))) / sum(pcomp(:,i)); % p is v_j - (x-a)/(yj-a)*v(x) for x=i'th target
      dv(j,i) = prefac(i) * ((s.x(j)-s.a)*p + (x(i)-s.x(j))*vb(j));
    end % pole-corrected for dv, gives bary stability for close node-targ pairs
  else  % interior:
    for l=1:numel(jj), j=jj(l); i=ii(l); % loop over node-targ pairs too close
      dv(j,i) = sum(comp(:,i).*(vb(j)-vb(:))) / sum(pcomp(:,i));
    end % bary for dv, gives bary stability for close node-targ pairs
  end
  vcp = prefac .* sum(dv.*comp.*invd) ./ J0; % bary form for deriv
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), j=jj(l); i=ii(l); % loop over hitting node-targ pairs
    vc(i) = vb(j);                     % replace each hit w/ corresp vb
    notj = [1:j-1, j+1:N];             % Schneider-Werner form for deriv @ node:
    vcp(i) = -sum(cw(notj).*(vb(j)-vb(notj))./(s.x(j)-s.x(notj)))/cw(j);
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vc vcp] = cauchycompmat_lsc2d(x,s,vb,side,o)
% Multiple vb-vector version of cauchycompeval_lsc2d, written by Gary Marple
% and Bowei Wu, 2014-2015. This can be called with eye(N) to fill the evaluation
% matrix (note that each column is not a holomorphic function, but by linearity
% when summed they give the right answer). When sent a single column vector
% for vb, this duplicates the action of cauchycompeval_lsc2d.
% Undocumented; notes by Barnett 6/12/16

% todo: check v' ext at the node (S-W form), seems to be wrong.

if nargin<5, o = []; end
if ~isfield(o,'delta'), o.delta = 1e-2; end
% dist param where $O(N^2.M) deriv bary form switched on.
% Roughly deriv errors are then limited to emach/mindist
% The only reason to decrease mindist is if lots of v close
% nodes on the curve with lots of close target points.

cw = s.cw;                    % complex speed weights
N = numel(s.x); M = numel(x);
n=size(vb,2);                  % # vb vectors

if nargout==1  % no deriv wanted... (note sum along 1-axis faster than 2-axis)
  comp = repmat(cw(:), [1 M]) ./ (repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1]));
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  %I0 = sum(repmat(vb(:),[1 M]).*comp); 
  I0 = permute(sum(repmat(permute(vb,[1 3 2]),[1 M 1]).*repmat(comp,[1 1 n]),1),[3 2 1]);
  J0 = sum(pcomp); % Ioakimidis notation
  vc = I0./(ones(n,1)*J0);                        % bary form
  if side=='e', vc = vc./(ones(n,1)*(x(:).'-s.a)); end         % correct w/ pole
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), vc(:,ii(l)) = vb(jj(l),:).'; end % replace each hit w/ corresp vb

else           % 1st deriv wanted...
  invd = 1./(repmat(s.x(:),[1 M]) - repmat(x(:).',[N 1])); % 1/displacement mat
  comp = repmat(cw(:), [1 M]) .* invd;
  if side=='e', pcomp = comp .* repmat(1./(s.x(:)-s.a), [1 M]);
  else pcomp = comp; end  % pcomp are weights and bary poles appearing in J0
  %I0 = sum(repmat(vb(:),[1 M]).*comp);
  I0 = permute(sum(repmat(permute(vb,[1 3 2]),[1 M 1]).*repmat(comp,[1 1 n]),1),[3 2 1]);
  J0 = sum(pcomp);
  if side=='e', prefac = 1./(x(:).'- s.a); else prefac = ones(1,M); end
  vc = (ones(n,1)*prefac) .* I0./(ones(n,1)*J0);                        % bary form (poss overall pole)
  %dv = repmat(vb(:),[1 M]) - repmat(vc(:).',[N 1]); % v value diff mat
    dv = repmat(permute(vb,[1 3 2]),[1 M 1]) - repmat(permute(vc,[3 2 1]),[N 1 1]); % v value diff mat
  [jj ii] = ind2sub(size(invd),find(abs(invd) > 1/o.delta)); % bad pairs indices
  if side=='e' % exterior:
    for l=1:numel(jj), j=jj(l); i=ii(l); % loop over node-targ pairs too close
      p = sum((comp(:,i)*ones(1,n)).*((ones(N,1)*vb(j,:))./(s.x(:)*ones(1,n)-s.a)-vb/(s.x(j)-s.a)),1) / sum(pcomp(:,i)); % p is v_j - (x-a)/(yj-a)*v(x) for x=i'th target
      dv(j,i,:) = prefac(i) * ((s.x(j)-s.a)*p + (x(i)-s.x(j))*vb(j,:));
    end % pole-corrected for dv, gives bary stability for close node-targ pairs
  else  % interior:
    for l=1:numel(jj), j=jj(l); i=ii(l); % loop over node-targ pairs too close
      dv(j,i,:) = sum((comp(:,i)*ones(1,n)).*(ones(size(vb,1),1)*vb(j,:)-vb),1) / sum(pcomp(:,i));
    end % bary for dv, gives bary stability for close node-targ pairs
  end
  % This is faster, but gives rounding errors when compared to the original
  % cauchycompeval.
  vcp = (ones(n,1)*(prefac./ J0)).* permute(sum(dv.*repmat(comp.*invd,[1 1 n]),1),[3 2 1]) ; % bary form for deriv
  
  % This vectorized form is slower, but does not give rounding errors.
  %vcp = (ones(n,1)*prefac).* permute(sum(dv.*repmat(comp,[1 1 n]).*repmat(invd,[1 1 n]),1),[3 2 1])./(ones(n,1)*J0) ; % bary form for deriv
  
  [jj ii] = ind2sub(size(comp),find(~isfinite(comp))); % node-targ coincidences
  for l=1:numel(jj), j=jj(l); i=ii(l); % loop over hitting node-targ pairs
    vc(:,i) = vb(j,:).';                     % replace each hit w/ corresp vb
    notj = [1:j-1, j+1:N];             % Schneider-Werner form for deriv @ node:
    vcp(i,:) = -sum((cw(notj)*ones(1,n)).*((ones(N-1,1)*vb(j,:))-vb(notj,:))./((s.x(j)-s.x(notj))*ones(1,n)),1)/cw(j);
  end
end
