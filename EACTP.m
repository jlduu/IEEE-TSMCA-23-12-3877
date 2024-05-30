% Function: Obtaining the complete fuzzy sociometric and DMs'weights
function Object=EACTP(x)
%x: the vector of trust values on a single path
d=size(x,2);  %the number of trust values in this path
% 1.Order effect factor
lam=zeros(d);
for k=1:d
    if d==1
       lam(k)=1;
    else 
    lam(k)=exp(-2.*(k-1)/((d-1).*d));
    end
% 2.OE-TP operator
p=[]; q=[];
if d==1
   Object=x;
else
    for k=1:d
        p=[p lam(k)*x(k)];
        q=[q (2-lam(k)*x(k))];
    end
    Object=2*prod(p)/(prod(q)+prod(p));
end
end
