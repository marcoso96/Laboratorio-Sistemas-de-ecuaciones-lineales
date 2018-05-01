function q=Buscapiv(fil, k)
  l=length(fil);
  q=k;
  grande=fil(k);
  for i=k:l
   
    if(abs(fil(i))>abs(grande)) 
      q=i;
      grande=fil(i);
    endif
  endfor
endfunction