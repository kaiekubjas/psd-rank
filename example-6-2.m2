restart
--psd rank
k = 3;
--matrix dimensions
p = 4;
q = 4;
R=QQ[a_(1,1,1)..a_(p,k,k),b_(1,1,1)..b_(q,k,k),x_(1)..x_(p*q)];

--ranks of matrices in the factorization
ranksA={1,1,1,1};
ranksB={1,1,1,1};

--construct A matrices
for i from 1 to p do (
    myList=for j from 1 to k list for l from 1 to ranksA#(i-1) list a_(i,j,l);
    A_(i)=matrix(myList)*transpose(matrix(myList));
)

--construct B matrices
for i from 1 to q do (
    myList=for j from 1 to k list for l from 1 to ranksB#(i-1) list b_(i,j,l);
    B_(i)=matrix(myList)*transpose(matrix(myList));
)
    
--construct list of entries of M
entriesList={flatten for i from 1 to p list for j from 1 to q list trace(A_(i)*B_(j))};
M=matrix(entriesList);
J=jacobian M;

--substitute random values for entries of A and B
substituteListA=flatten flatten for i from 1 to p list for j from 1 to k list for l from 1 to k list a_(i,j,l)=>random(1,100);
substituteListB=flatten flatten for i from 1 to q list for j from 1 to k list for l from 1 to k list b_(i,j,l)=>random(1,100);
substituteList=join(substituteListA, substituteListB);
J2=sub(J,substituteList);
rank(J2)
