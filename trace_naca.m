load solutionnaca.dat;
m=load_gmsh('wing_naca.msh');
clf;
xP=zeros(3,m.nbTriangles);
yP=xP;
zP=xP;
for j=1:m.nbTriangles
    for k=1:3
        noe=m.TRIANGLES(j,k);
        xP(k,j)=m.POS(noe,1);
        yP(k,j)=m.POS(noe,2);
        zP(k,j)=solutionnaca(noe);
    end
end
valtrix=sum(xP)'/3;
valtriy=sum(yP)'/3;
valtriu=sum(zP)'/3;
pp=[m.POS(1:m.nbNod,1)';m.POS(1:m.nbNod,2)'];
uu=valtriu;
tt=m.TRIANGLES(1:m.nbTriangles,:)';
valneu=pdeprtni(pp,tt,uu');
figure(1);
pdecont(pp,tt,valneu,100);

axis tight;

