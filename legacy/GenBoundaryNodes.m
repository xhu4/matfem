function boundary_nodes = GenBoundaryNodes(nx, ny)

% nbn = (nx+ny)*2-4;
idx_map = reshape(1:nx*ny,[ny,nx]);

bottom = idx_map(1,2:end);
nbottom = numel(bottom);

right = idx_map(1:end, end)';
nright = numel(right);

top = flip(idx_map(end,1:end-1));
ntop = numel(top);

left = flip(idx_map(1:end-1, 1))';
nleft = numel(left);

boundary_nodes = [1*ones(1,nbottom), 1*ones(1,nright), 1*ones(1,ntop), 1*ones(1,nleft);
	bottom, right, top, left];