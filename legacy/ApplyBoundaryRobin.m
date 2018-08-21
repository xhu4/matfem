function [R, w] = ApplyBoundaryRobin(be, Pb, Tb, nqpts, dxdy_R, dxdy_w, ...
	f_R, f_w, basis_type)

Nlb = size(Tb, 1);
Nb = size(Pb,2);

w = zeros(Nb, 1);
R = zeros(Nb, Nb);

for iEdge = find(be(1,:)==3)
	iElem = be(2, iEdge);
	nodes = Pb(:,be(3:4, iEdge))';
	idx = Tb(:, iElem);
	vElem = Pb(:,idx(1:3))';
	[x, weights] = lgwt(nqpts, nodes(1,1), nodes(2,1));
	[y, ~] = lgwt(nqpts, nodes(1,2), nodes(2,2));
	
	for beta = 1:Nlb
		if ~isempty(f_w)
			r = IntegrateGauss2d(vElem, beta, dxdy_w, f_w, x, y, weights, 1, basis_type);
			w(idx(beta),1) = w(idx(beta),1) + r;
		else
			w = [];
		end
		
		if ~isempty(f_R)
			for alpha = 1:Nlb
				r = IntegrateGauss2d(vElem, [beta, alpha], dxdy_R, f_R, x, y, weights, 1, basis_type);
				R(idx(beta), idx(alpha)) = R(idx(beta), idx(alpha)) + r;
			end
		else
			R = [];
		end
	end
end