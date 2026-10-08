function p=oosetfemops(p)
    gr=p.pdeo.grid;
    [K,M,~]=p.pdeo.fem.assema(gr,1,1,1); 
    p.mat.M=kron([[1,0];[0,1]],M); 
    p.mat.K=K;
end