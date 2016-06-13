function h = showsegment(s, trlist)
% SHOWSEGMENT   plot segment(s) & possibly translated copies
%
% h = showsegment(s) where s is segment struct with s.x nodes, optionally s.nx
%  normal vectors, adds to the current axes a plot of the segment.
%  If s is a cell array of segment structs, it plots all of them.
%
% h = showsegment(s, trlist) also plots translated copies given by list of
%  complex numbers in trlist.
%
% No arguments does a self-test.

% Based on variants of showseg.  Barnett 6/12/16

if nargin<1, test_showsegment; return; end
if nargin<2, trlist = 0; end
if iscell(s), for i=1:numel(s), showsegment(s{i},trlist); end, return, end
hold on
for i=1:numel(trlist)
  plot(s.x+trlist(i), 'b.-');
  if isfield(s,'nx')
    l=0.05; plot([s.x, s.x+l*s.nx].'+trlist(i), 'k-');
  end
end
axis equal xy

%%%%%%%
function test_showsegment
N = 100; s.x = exp(2i*pi*(1:N)/N);
figure; showsegment(s);            % plain, no normals
s = setupquad(s);                  % adds normals
[xx yy] = meshgrid([-3 0 3]);
trlist = xx(:)+1i*yy(:);
figure; showsegment(s,trlist);      % check translation copies
s2 = s; s2.x = s.x*0.5;            % smaller circles
figure; showsegment({s s2},trlist);  % check cell array
