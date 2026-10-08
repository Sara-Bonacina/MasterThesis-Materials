function out=ToxicityX1Dbra(p,u)
    out=[u(p.np+1);... % 2nd component in 0
         p.Om*sum(p.mat.M(1:p.np,1:p.np)*abs(u(1:p.np)));...% |u|_L^1
         p.Om*sum(p.mat.M(1:p.np,1:p.np)*abs(u(p.np+1:2*p.np)))];% |v|_L^1
end
