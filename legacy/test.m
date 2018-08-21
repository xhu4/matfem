import MatFem.*

nx = 3;
ny = 3;
mesh = RectMesh([0,1,0,1],[nx, ny]);
assert(all(mesh.x == linspace(0,1,nx)));
assert(mesh.ny == ny);
assert(mesh.hy == 1/(ny-1));

spc = RectFcnSpc(mesh, 2);
Tb = [ 1  3  3  5 11 13 13 15;
	  11 11 13 13 21 21 23 23;
	   3 13  5 15 13 23 15 25;
	   6  7  8  9 16 17 18 19;
	   7 12  9 14 17 22 19 24;
	   2  8  4 10 12 18 14 20];
x = 0:.25:1;
Pb = [repelem(x, 5);
	  repmat(x, 1, 5)];
assert(isequal(spc.Tb, Tb));
assert(isequal(spc.Pb, Pb));