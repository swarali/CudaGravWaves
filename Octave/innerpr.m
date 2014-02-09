function y= innerpr(h1,h2,fr_l,fr_h)

#df1 = (fr_h -fr_l)/size(h1)(2);
df = 1/size(h1)(2);
y1 = conj(h1) * conj(h2)'
y2 = conj(h2) * conj(h1)'
#imag(h2)
#imagine=sum(imag(h2).^2)/sum(real(h2).^2)
#imagine=sum(imag(h2).^2)


y = (y1+y2)*df/2;


endfunction
