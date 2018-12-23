function [p, f] = dps1()
p.dim = 2;
p.eta = 1; p.alpha = 1;
p.phi_m = 1; p.phi_f = 1;
p.k_m = .01; p.k_f = 1;
p.mu = 1; p.nu = 1; p.rho = 1; p.sigma = 1;
p.C_mt = 1; p.C_ft = 1;
p.PI = p.k_f * eye(p.dim);

f.pm = @(x,y,t) sin(x.*y.^2-y.^3).*cos(t);
f.lap_pm = @(x,y,t) cos(t).*(2*(x-3*y).*cos(y.^2.*(x-y))- ...
	2*y.^2.*(2*x.^2-6*x.*y+5*y.^2).*sin(y.^2.*(x-y)));
f.dpmdt = @(x,y,t) -sin(x.*y.^2-y.^3).*sin(t);
f.dpmdx = @(x,y,t) cos(x.*y.^2-y.^3).*y.^2.*cos(t);
f.dpmdy = @(x,y,t) cos(x.*y.^2-y.^3).*(2*x.*y-3*y.^2).*cos(t);

f.pf = @(x,y,t) (2-pi*sin(pi*x)).*(cos(pi*(1-y))-y).*cos(2*pi*t);
f.dpfdt = @(x,y,t) -2*pi*(2-pi*sin(pi*x)).*(cos(pi*(1-y))-y).*sin(2*pi*t);
f.dpfdx = @(x,y,t) pi^2*cos(2*pi*t).*cos(pi*x).*(y+cos(pi*y));
f.dpfdy = @(x,y,t) cos(2*pi*t).*(2-pi*sin(pi*x)).*(pi*sin(pi*y)-1);
f.lap_pf = @(x,y,t) pi^2*cos(2*pi*t).*((2-2*pi*sin(pi*x)).*cos(pi*y)-pi*y.*sin(pi*x));

f.uc1 = @(x,y,t) (x.^2.*y.^2+exp(-y)).*cos(2*pi*t);
f.duc1dt = @(x,y,t) -2*pi*(x.^2.*y.^2+exp(-y)).*sin(2*pi*t);
f.duc1dx = @(x,y,t) (2*x.*y.^2).*cos(2*pi*t);
f.duc1dy = @(x,y,t) (2*x.^2.*y-exp(-y)).*cos(2*pi*t);
f.duc1dxx = @(x,y,t) (2*y.^2).*cos(2*pi*t);
f.duc1dxy = @(x,y,t) 4*x.*y.*cos(2*pi*t);
f.duc1dyy = @(x,y,t) (2*x.^2+exp(-y)).*cos(2*pi*t);

f.uc2 = @(x,y,t) (-2/3*x.*y.^3+(2-pi*sin(pi*x))).*cos(2*pi*t);
f.duc2dt = @(x,y,t) -2*pi*(-2/3*x.*y.^3+(2-pi*sin(pi*x))).*sin(2*pi*t);
f.duc2dx = @(x,y,t) (-2/3*y.^3-pi^2*cos(pi*x)).*cos(2*pi*t);
f.duc2dy = @(x,y,t) (-2*x.*y.^2).*cos(2*pi*t);
f.duc2dyy = @(x,y,t) -4*x.*y.*cos(2*pi*t);
f.duc2dxy = @(x,y,t) -2*y.^2.*cos(2*pi*t);
f.duc2dxx = @(x,y,t) pi^3*sin(pi*x).*cos(2*pi*t);

f.p = @(x,y,t) (pi*sin(pi*x)-2).*cos(2*pi*y).*cos(2*pi*t);
f.dpdx = @(x,y,t) pi^2*cos(pi*x).*cos(2*pi*y).*cos(2*pi*t);
f.dpdy = @(x,y,t) -2*pi*(pi*sin(pi*x)-2).*sin(2*pi*y).*cos(2*pi*t);

f.Q = @(x,y,t) p.sigma*p.k_m/p.mu*(f.pm(x,y,t)-f.pf(x,y,t));
f.Qplus = @(x,y,t) p.phi_m*p.C_mt*f.dpmdt(x,y,t)-p.k_m/p.mu*f.lap_pm(x,y,t)+f.Q(x,y,t);

f.f1 = @(x,y,t) -2*pi*(x.^2.*y.^2+exp(-y)).*sin(2*pi*t)- ...
	p.nu*(2*y.^2+2*x.^2+exp(-y)).*cos(2*pi*t)+pi^2*cos(pi*x).*cos(2*pi*y).*cos(2*pi*t);
f.f2 = @(x,y,t) -2*pi*(-2/3*x.*y.^3+2-pi*sin(pi*x)).*sin(2*pi*t)- ...
	p.nu*cos(2*pi*t).*(pi^3*sin(pi*x)-4*x.*y)-2*pi*(pi*sin(pi*x)-2).*sin(2*pi*y).*cos(2*pi*t);

f.qp = @(x,y,t) p.phi_f*p.C_ft*f.dpfdt(x,y,t)-p.k_f/p.mu*f.lap_pf(x,y,t)-f.Q(x,y,t);

testfcns(p, f);
end

function testfcns(p, f)
% on \Omega_D (dual porosity subdomain)
[x,y,t] = meshgrid(0:.1:1,0:.05:0.75,0:.1:1);
% 2.1
lhs = p.phi_m*p.C_mt*f.dpmdt(x,y,t)-p.k_m/p.mu*f.lap_pm(x,y,t);
rhs = -f.Q(x,y,t)+f.Qplus(x,y,t);
almosteq(lhs, rhs);

% 2.2
lhs = p.phi_f*p.C_ft*f.dpfdt(x,y,t)-p.k_f/p.mu*f.lap_pf(x,y,t);
rhs = f.Q(x,y,t)+f.qp(x,y,t);
almosteq(lhs, rhs);

% on \Omega_C (conduit subdomain)
[x,y,t] = meshgrid(0:.1:1,-.25:.05:0,0:.1:1);
% 2.3
lhs = f.duc1dt(x,y,t)-(2*p.nu*f.duc1dxx(x,y,t)-f.dpdx(x,y,t)+...
	p.nu*(f.duc1dyy(x,y,t)+f.duc2dxy(x,y,t)));
rhs = f.f1(x,y,t);
almosteq(lhs, rhs);

lhs = f.duc2dt(x,y,t)-(2*p.nu*f.duc2dyy(x,y,t)-f.dpdy(x,y,t)+...
	p.nu*(f.duc1dxy(x,y,t)+f.duc2dxx(x,y,t)));
rhs = f.f2(x,y,t);
almosteq(lhs, rhs);

% 2.4
lhs = f.duc1dx(x,y,t)+f.duc2dy(x,y,t);
almosteq(lhs,0);

% on \Gamma_{CD} (interface)
[x,y,t] = meshgrid(0:.1:1,0,0:.1:1);
% 2.5
lhs = f.dpmdy(x,y,t);
almosteq(lhs,0);

% 2.6
lhs = f.uc2(x,y,t);
rhs = -p.k_f/p.mu*f.dpfdy(x,y,t);
almosteq(lhs,rhs);

% 2.7
lhs = -2*p.nu*f.duc2dy(x,y,t)+f.p(x,y,t);
rhs = f.pf(x,y,t)/p.rho;
almosteq(lhs,rhs);

% 2.8
lhs = -p.nu*(f.duc1dy(x,y,t)+f.duc2dx(x,y,t));
rhs = p.alpha*p.nu*sqrt(p.dim)/sqrt(trace(p.PI))*(f.uc1(x,y,t)+p.k_f/p.mu*f.dpfdx(x,y,t));
almosteq(lhs,rhs);

lhs = -p.nu*2*f.duc2dy(x,y,t);
rhs = p.alpha*p.nu*sqrt(p.dim)/sqrt(trace(p.PI))*(f.uc2(x,y,t)+p.k_f/p.mu*f.dpfdy(x,y,t));
almosteq(lhs,rhs);
end


function almosteq(a,b,tol)
if nargin < 3
	tol = 1e-8;
end
d = abs(a-b);
assert(all(d(:)<tol));
end