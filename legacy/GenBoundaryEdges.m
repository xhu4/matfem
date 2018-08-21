function boundaryEdges = GenBoundaryEdges(range, h_x, h_y, T)
% GENBOUNDARYEDGES generates index matrix for boundary edges
%
%	See also GENBOUNDARYNODES

left = range(1);
right = range(2);
bottom = range(3);
top = range(4);

if nargin<4
	[~,T] = Gen2dTriMesh([left, right ,bottom, top], h_x, h_y);
end

assert(mod(right-left, h_x)==0 && mod(top-bottom, h_y)==0)

nx = (right-left)/h_x;
ny = (top-bottom)/h_y;
n_ele = size(T,2);

% nbe = 2*(nx+ny);

bottom_eles = 1:ny*2:n_ele;
bottom_idx = T([1,2],bottom_eles);
nbottom = numel(bottom_eles);

right_eles = (nx-1)*ny*2+2:2:n_ele;
right_idx = T([2,3],right_eles);
nright = numel(right_eles);

top_eles = flip(bottom_eles)+1+ny;
top_idx = T([3,1],top_eles);
ntop = numel(top_eles);

left_eles = ny*2-1:-2:1;
left_idx = T([3,1],left_eles);
nleft = numel(left_eles);

boundaryEdges = [1*ones(1,nbottom), 1*ones(1,nright), 1*ones(1,ntop), 1*ones(1,nleft);
				  bottom_eles, right_eles, top_eles, left_eles;
				  bottom_idx, right_idx, top_idx, left_idx];