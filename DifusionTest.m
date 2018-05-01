# Test Ejercicio 2
# Marcos Obando

tol=1e-6;

for n=1:14
  clear L U A;
  N(n)=8*n^3;  
  [A, b]=Difusion(n);
  
  #metodo LU
  [L,U]=lu(A); 
  
  t1=time();
  
  #comparar con y sin tomar en cuenta la descomposicion
  y=b\L;
  x=y'\U;
  t2_LU(n)=time()-t1;
  
  clear L U
  #Metodo sin condicionar pcg
  t1=time();
  
  [x,~,~, iter(n)]=pcg(A, b,tol, N(n));
  
  t2(n)=time()-t1;
  
  #metodo pcg condicionado
  [L, U]=ilu(A);
  
  t1=time();
    
  [x,~,~,iter_cond(n)]=pcg(A, b, tol,MAXIT=8*n^3,L,U);
  
  t2_Cond(n)=time()-t1;
  
  #Metodo barra\
  
  t1=time();
  
  x=A\b;
  t2_barra(n)=time()-t1;
  
endfor

for n=1:6

  N_it(n)=8*n^3;  
  [A, b]=Difusion(n);
  
  cond_A(n)=cond(A);
  [L,U]=ilu(A);
  
  C_inv=inv(sqrtm(L*U));
  cond_C(n)=cond(C_inv*A*C_inv);
end

clear N;

for n=1:14
    
  N(n)=8*n^3; 

end


figure(1);
loglog(N, t2, 'ko', 'markersize', 13,'markerfacecolor', 'b',  N,t2_Cond,'kdiamond','markersize', 13,'markerfacecolor', 'y',N, t2_LU, 'ko','markersize', 13,'markerfacecolor', 'r',N, t2_barra,'ko', 'markersize', 13, 'markerfacecolor', 'k');
h=legend('Grad. Conj. sin precondicion','Grad. Conj. precondicionado' , 'Factorizacion LU', 'Barra Invertida','location', 'northwest');
set(h,'fontsize', 18);
set(gca(), 'fontsize', 26);
xlabel('Dimension de la matriz', 'fontname', 'Helvetica', 'fontsize', 28);
ylabel('Tiempo de calculo[s]', 'fontname', 'Helvetica', 'fontsize', 28);


figure(2);
loglog(N_it,cond_A, 'ro','linestyle', 'none', 'markersize', 13, 'markerfacecolor', 'b', 
N_it, cond_C, 'bsquare','linestyle', 'none', 'markersize', 13, 'markerfacecolor', 'r');
h=legend('Sistema sin condicionar','Sistema precondicionado', 'location', 'northwest');
set(h,'fontsize', 26);
set(gca(), 'fontsize', 26);
xlabel('Dimension del sistema a resolver', 'fontname', 'Helvetica', 'fontsize', 28);
ylabel('Numero de condicion del sistema', 'fontname', 'Helvetica', 'fontsize', 28);
