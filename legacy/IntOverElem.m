function ret = IntOverElem(vElem, f, nqpts)
% INTOVERELEM integrates function f over element vElem
%
%	See also INTOVERMESH

[X, Y, Wx, Wy] = triquad(nqpts, vElem);
ret = Wx'*f(X,Y)*Wy;