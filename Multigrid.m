clear;
clc;
lev=6;
epsilon=1;
alfa=1;
Nmax=100000;
n1=2;
n2=1;
omega=6/7;
mu=1;
tol=1e-10;


for i=1:lev
    n=2^(lev-i+1);
    T{i}=gallery('tridiag',n-1)*n*n;
    I{i}=speye(n-1,n-1);
    C{i}=kron(I{i},T{i})+kron(T{i},I{i});
    K{i}=sparse(n*diag(-1*ones((n-1)^2,1)) + n*diag(ones((n-1)^2-1,1),-1));
    b{i}=zeros((n-1)^2,1);
    A{i}=-epsilon*C{i}+alfa*K{i};  
    D{i}=diag(A{i});
    x{i}=zeros((n-1)^2,1);
end

for i=1:lev-1
    nf=2^(lev-i+1);
    nc=2^(lev-i);
    PP{i}=spalloc(nf-1,nc-1,2*(nf-1));
    for j=1:nc-1
        PP{i}(2*(j-1)+1:2*(j)+1,j)=[.5 1 .5]';
    end
    P{i}=kron(PP{i},PP{i});
    R{i}=0.25*P{i}';
end

b{1}=ones((2^lev-1)^2,1);
normb=norm(b{1});

tic
for i=1:Nmax
    x=W(A,D,P,R,x,b,n1,n2,omega,1,lev,mu);% W-cycle 
    rc=norm(b{1}-A{1}*x{1});
     disp(rc);
    if rc<tol*normb
        disp(['Multigrid mu(' num2str(mu) ')-Cycle converged at ' num2str(i) ' iterations']);
        break;
    end
end
toc
