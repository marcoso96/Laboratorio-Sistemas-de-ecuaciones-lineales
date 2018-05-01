t_fin=zeros(100);

for i = 1:100

  A=rand(i,i);
  x=rand(i,1);

  t=time();
  [L, U]=decomp_lu(A);
  t_fin(i)=time()-t;
  
end 

figure(1)
loglog(t_fin, 'bo','linestyle', 'none', 'markersize', 16, 'markerfacecolor', 'b');
set(gca(), 'fontsize', 26);
xlabel('Dimension de la matriz A', 'fontsize', 30,'fontname', 'Helvetica');
ylabel('Tiempo de calculo [s]','fontsize', 30, 'fontname', 'Helvetica');
