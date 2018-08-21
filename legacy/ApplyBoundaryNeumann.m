function v = ApplyBoundaryNeumann(be, Pb, Tb, nqpts, dxdy, f, basis_type)

Nlb = size(Tb, 1);
Nb = size(Pb,2);

v = zeros(Nb, 1);

for iEdge = find(be(1,:)==2)
	iElem = be(2, iEdge);
	nodes = Pb(:,be(3:4, iEdge))';
	vElem = Pb(:,Tb(1:3,iElem))';
	[x, w] = lgwt(nqpts, nodes(1,1), nodes(2,1));
	[y, ~] = lgwt(nqpts, nodes(1,2), nodes(2,2));
	
	for beta = 1:Nlb
		r = IntegrateGauss2d(vElem, beta, dxdy, f, x, y, w, 1, basis_type);
		v(Tb(beta, iElem),1) = v(Tb(beta, iElem),1) + r;
	end
end