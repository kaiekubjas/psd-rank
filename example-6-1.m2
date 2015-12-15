restart
r = 3;
k = 4;
p = 4;
q = 4;
R=QQ[x,y,X,Y,a_(1)..a_(r)];

--circulant matrix
C=matrix{{x, y, 1, y},{ y ,x, y, 1}, {1, y, x, y}, {y, 1, y, x}};
C2=matrix{(flatten entries C)/(x->x^2)}

--construct A matrices
A_(1)=matrix{{1,0,0},{0,0,0},{0,0,0}};
A_(2)=matrix{{0,0,0},{0,1,0},{0,0,0}}
A_(3)=matrix{{0,0,0},{0,0,0},{0,0,1}}
for i from p to p do (
    myList=for j from 1 to r list {a_(j)};
    A_(i)=matrix(myList)*transpose(matrix(myList));
)

--construct B matrices
for i to q-1 do (
    myList=for j to r-1 list for l to r-1 list C_(j,i)*C_(l,i);
    B_(i+1)=matrix(myList)
)
    
--construct list of entries of M
entriesList={flatten for i from 1 to p list for j from 1 to q list trace(A_(i)*B_(j))};
M=matrix(entriesList);

--construct the ideal of the circulant matrix and eliminate variables of the factorization
myList=for j to 15 list C2_(0,j)-M_(0,j);
I=ideal(myList)+ideal(X-x^2,Y-y^2)
J=eliminate({a_(1),a_(2),a_(3)},I)
K=eliminate({x,y},J)
toString factor (entries gens  K)#0#0
