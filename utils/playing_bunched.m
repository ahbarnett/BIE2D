% play w/ reparam circle to handle singularities. Barnett 9/27/17
clear
N=1e3;
t=(1:N)/N * 2*pi; h = 2*pi/N;
be = 10;
yp = cosh(be*sin(2*t));
I = sum(yp)*h;
al = 2*pi/I;
yp = yp*al;
figure; plot(t,yp,'.-'); title('deriv y''(t)')
y = perispecint(yp);
figure; plot(t,y,'.-'); title('reparam map y(t)')
figure; plot(y,0*y,'.'); title('quadr nodes on [0,2pi)')

r = 0.5; s = wobblycurve(r,0,1,N);
s.x = s.Z(y);
s.xp = yp.*s.Zp(y);
s=rmfield(s,{'Z','Zp','Zpp'});  % otherwise get used
s = setupquad(s);  % generate from s.x and s.xp
s.cur = 2*ones(N,1);     % for circle only
figure; showsegment(s)


