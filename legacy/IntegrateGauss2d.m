function result = IntegrateGauss2d(vE, bases, dxdy, c, X, Y, Wx, Wy, basis_type)
% INTEGRATEGAUSS2D integrates using gauss quadrature method.
%	INTEGRATEGAUSS2D evaluates:
%		\int_E F dxdy
%	E is the element with vertices vE.
%	F is the product of derivatives of all bases function, derivatives of
%	each base function is determined by each row of dxdy.
%
% See also ASSEMBLESTIFF2D, ASSEMBLELOAD2D

assert(all(numel(bases) == size(dxdy, 1)));

result = Wx*Wy';

for i = 1:length(bases)
	r = Eval2dTriBasis(X, Y, vE, bases(i), basis_type, dxdy(i,:));
	result = result .* r;
end

coeff_x = c(X,Y);

result = sum(sum((result.*coeff_x)));