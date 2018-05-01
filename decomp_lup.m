function [l, u,p] = decomp_lup(a)
  n=length(a);
  l=zeros(n,n);
  p=eye(n,n);
  u=a;
  
  for k=1:(n-1)#columnas que escalona
    
    q=Buscapiv(u(:,k), k); #busco el mayor pivote en la columna k-esima
    
    temp=u(k,:);
    u(k,:)=u(q,:);
    u(q,:)=temp;
    
    temp=l(k,:);
    l(k,:)=l(q,:);
    l(q,:)=temp;
    
    temp=p(k,:);
    p(k,:)=p(q,:);
    p(q,:)=temp;
       
    for i=(k+1):n #elementos bajo la diagonal
      
      l(i,k)=u(i,k)/u(k,k);
      
      for j=k:n
        
        u(i,j)=u(i,j)-l(i,k)*u(k, j);
        
      endfor
    endfor
  endfor 
  l=eye(n,n)+l;
endfunction  