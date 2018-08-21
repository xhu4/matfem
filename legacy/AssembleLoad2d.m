function b = AssembleLoad2d(Pb, Tb, nqpts, dxdy, f, basis_type)
% ASSEMBLELOAD2D assembles the load vector.
%	ASSEMBLELOAD2D builds the vector b with
%		b_i = \int f phi_i?
%	with derivative of phi_i determined by dxdy.
%
% See also ASSEMBLESTIFF2D, INTEGRATEGAUSS2D
[Nlb, Ne] = size(Tb);
Nb = size(Pb,2);

b = zeros(Nb, 1);

for iE = 1:Ne
	vE = Pb(:,Tb(1:3,iE))';
	[X, Y, Wx, Wy] = triquad(nqpts, vE);
	
	for beta = 1:Nlb
		r = IntegrateGauss2d(vE, beta, dxdy, f, X, Y, Wx, Wy, basis_type);
		b(Tb(beta, iE),1) = b(Tb(beta, iE),1) + r;
	end
end