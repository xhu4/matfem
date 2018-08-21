function r = intOverElem(spc, iE, ders, coefFcn, nqpts)
vE = spc.elem(iE);
vE = vE(1:3,:);
r = MatFem.intOver(vE, spc.basisType, ders, coefFcn, nqpts);
end