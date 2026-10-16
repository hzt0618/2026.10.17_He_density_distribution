function res=Pnm(n,m,sita)

if m==0 && n==0
    res=1;
end

if m==0 && n==1
    res=cos(sita);
end

if m==1 && n==0
    res=0;
end

if m==1 && n==1
    res=sin(sita);
end

if m==n && n>1
    res=(2*n-1)*sin(sita).*Pnm(n-1,m-1,sita)/sqrt((n+m)*(n+m-1));
end

if m~=n && n>1
    res=((2*n-1)*cos(sita).*Pnm(n-1,m,sita)-sqrt((n-m-1)*(n+m-1))* ...
        Pnm(n-2,m,sita))/sqrt((n+m)*(n-m));
end

if m>n
    res=0;
end

end

