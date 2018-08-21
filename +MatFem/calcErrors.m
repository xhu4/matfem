function T = calcErrors(fcnSpc, u, utrue, uders, h)
% CALCERRORS calculate the errors for our fem solution
%
% CALCERRORS(fcnSpc, u, utrue, uders, h)

if nargin == 0
	T = table({},[],[],[]);
	T.Properties.VariableNames={'h', 'error_inf', 'error_L2', 'error_H1'};
	return
end

hx = h;

error_L2 = MatFem.errorNorm(fcnSpc, u, utrue, 'L2');

error_inf = MatFem.errorNorm(fcnSpc, u, utrue, 'inf');

error_H1 = MatFem.errorNorm(fcnSpc, u, uders, 'H1');

h = {strtrim(rats(hx))};
T = table(h, error_inf, error_L2, error_H1);
