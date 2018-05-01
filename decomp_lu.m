function [l, u] = decomp_lu(a)
  n=length(a);
  l=eye(n,n);

  for k=1:(n-1)#columnas que escalona
    
    for i=(k+1):n #elementos bajo la diagonal
      
      l(i,k)=a(i,k)/a(k,k);
      
      for j=k:n
      
        a(i,j)=a(i,j)-l(i,k)*a(k, j);
      
      endfor
    endfor
  endfor

    u=a;
endfunction  
