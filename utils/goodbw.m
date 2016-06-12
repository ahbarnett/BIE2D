function c = goodbw(n)
% Returns a colormap of n-by-3 size, that reproduces well in B/W.
% If n>0, n colors increase from dark to light; if n<0, then (-n) colors
% go other way.
% adapted from Carey Rappaport (C) 2002, CMR documentation.
% Barnett 9/26/12

%CMRmap=[0 0 0;.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1]; % the original
%CMRmap = CMRmap .* repmat([1 1.5 1.45 1.4 1.35 1.3 1.25 1.2 1], [3 1])';
%CMPmap = min(max(CMRmap,0),1) % clip to [0,1]

%CMRmap=[0 0 0;.2 .2 .7;.5 .2 .9;.8 .3 .6;1 .5 .2;1 .7 .4;1 .9 .4;1 1 .7;1 1 1];  % brightened up a bit!  good linearity, red is bit weak

CMRmap=[0 0 0;.2 .2 .7;.5 .2 .9;.8 .3 .6;1 .45 .3;1 .7 .3;1 .9 .4;.9 1 .8;1 1 1]; % brightened up a bit!  good linearity, better, a touch of green at the bottom

if nargin<1, n = -256; end  % default

x = 1:8/(abs(n)-1):9;	% n color levels instead of 9
x1 = 1:9;
for i = 1:3
  c(:,i) = spline(x1,CMRmap(:,i),x)'; % spline fit intermediate values
end
%c = min(max(c,0),1); % clip to [0,1]
c = abs(c/max(max(c)));           % eliminate spurious values outside of range
c(end,:) = [1 1 1];   % ensure totally white to start
if n<0, c = c(end:-1:1,:); end   % reverse order
