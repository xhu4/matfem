function ret = EvalFcn(u, X, Y, Pb, Tb, iE, basisType, diff)
% EVALFCN evaluates function u = sum u_i phi_i at point(s) (X, Y)

Nlb = size(Tb, 1);
vE = VtcsOfElem(Pb, Tb, iE);

ret = 0;
for ib = 1:Nlb
	iv = Tb(ib, iE);
	ret = ret + u(iv)*Eval2dTriBasis(X, Y, vE, ib, basisType, diff);
end