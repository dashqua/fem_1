lc=0.05; //importance physique du noeud
nw=20; //nb de points sur chaque tranche d'aile

// Definition des points
Point(1) = {-1,-1,0,lc};
Point(2) = {2,-1,0,lc};
Point(3) = {2,1,0,lc};
Point(4) = {-1,1,0,lc};

// Points Wing
For i In {1:nw}
	t=(i-1)/(nw-1);
	Point(4+i) = {t,0.17735*Sqrt(t)-0.075597*t-0.212836*t*t+0.17363*t*t*t-0.06254*t*t*t*t,0,lc};
EndFor
For i In {1:nw-3}
	t=(nw-i-2)/(nw-2);
	Point(4+nw+i) = {t,-(0.17735*Sqrt(t)-0.075597*t-0.212836*t*t+0.17363*t*t*t-0.06254*t*t*t*t),0,lc};
EndFor

// Definition des lignes
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3 ,4};
Line(4) = {4, 1};

// LinesWing
For i In {1:2*nw-4}
	Line(4+i)={4+i,4+i+1};
EndFor
Line(2*nw+1)={2*nw+1,5};

// Definition du domaine
Line Loop(1)={1,2,3,4};
Line Loop(2) = {5:1+nw*2};
Plane Surface(1) = {1,2};

// Definitions des zones physiques avec attribution d'un identifiant
Physical Line(1001) = {1, 2, 3, 4};
Physical Line(1002) = {5:1+nw*2};
Physical Surface(1003) = {1};

