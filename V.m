function x=V(A,D,P,R,x,b,n1,n2,omega,l,lmax,mu)
if l==lmax
    x{l}=A{l}\b{l};
else
    L{l} = tril(A{l});
    TGS{l}=-L{l}\triu(A{l},1);
    for i=1:n1
%         x{l}=x{l}+omega*(b{l}-A{l}*x{l})./D{l};
        x{l}=TGS{l}*x{l}+L{l}\b{l};
    end
    b{l+1}=R{l}*(b{l}-A{l}*x{l});
    x{l+1}=zeros(length(x{l+1}),1);
    for i=1:mu
        x=mgmu(A,D,P,R,x,b,n1,n2,omega,l+1,lmax,mu);
    end
    x{l}=x{l}+P{l}*x{l+1};
    for i=1:n1
%         x{l}=x{l}+omega*(b{l}-A{l}*x{l})./D{l};
         x{l}=TGS{l}*x{l}+L{l}\b{l};
    end
end
end