function vE = VtcsOfElem(Pb, Tb, iE)
% GETVTCSOFELEM get the coordinates of all vertices of the iE-th element

vE = Pb(:,Tb(1:3,iE))';