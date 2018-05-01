for i = 1:80

  A=rand(i,i);
  x=rand(i,1);
  A(1:1)=1e-10;
   
  b=A*x;
  
  [L, U, P_s]=decomp_lup(A);
  [L_s, U_s]=decomp_lu(A);
  
  y_lu=L_s\b;
  X_lu=U_s\y_lu;
  
  Y_lup=L\(P_s*b);
  X_lup=U\Y_lup;
  
  err_LUP(i)=norm(x-X_lup);
  err_LU(i)=norm(x-X_lu);
end

for i = 1:80  
  prueba(i)=i;
end 

figure(1);
loglog(prueba, err_LU, 'ko', 'markersize', 16,'markerfacecolor', 'b',  prueba,err_LUP,'kdiamond','markersize', 16,'markerfacecolor', 'y');
h=legend('LU','LUP' ,'location', 'northwest');
set(h,'fontsize', 26);
set(gca(), 'fontsize', 26);
xlabel('Dimension de la matriz', 'fontname', 'Helvetica', 'fontsize', 30);
ylabel('Tiempo de calculo[s]', 'fontname', 'Helvetica', 'fontsize', 30);