function chirpfr5(fr,n,samp,fr_end,inc)
	#plot the correlation values of template and signal at various phases
	#fig1 - chirp sig in time domain		fig2 - ftt of chirp signal
	#{ 
	fr - template fr, 
	n - duration of the signal, 
	samp - sampling frequency, 
	fr_end - end frequency,
	k - increase in fr per sample of chirp signal, k= (fr_end-fr)/(n*samp+1)
	t_n - duration of the template
	#}

	k= (fr_end-fr)/(n*samp+1);
	x=linspace(0,n, (n*samp)+1);
	#generating chirp signal 
	sig= zeros(1,n*samp+1);
	fr_1=(x.*k)+fr;
	sig(1:n*samp +1) = x.*sin(x.*(fr_1.*2*pi));
	plot(x,sig,"-",x,x.*0,"-");	
	print -dsvg fig1.svg;

	sig_nor = sig./sqrt(sum(sig.^2));
	sigfr = fft(sig_nor);
	sigfrpl =abs(sigfr(1:n*samp/2));
	plot(x(1:n*samp/2).*samp/n,sigfrpl,"-");
	print -dsvg fig2.svg;

	max_cor = innerpr(sigfr,sigfr,0,samp/2)

	#changed frequency
	temp_inc = 0;
	tempcor=1.00;

	while(tempcor>0.97)
		temp_inc = temp_inc - inc;
		sig_new(1:n*samp +1) = x.*sin(x.*((fr_1+ temp_inc).*2*pi));
	
	#	plot(x,sig_new,"-",x,x.*0,"-");
	#	print -dsvg fig3.svg;
	
		signew_nor = sig_new./sqrt(sum(sig_new.^2));
		datafr = fft(signew_nor);
		datafrpl =abs(datafr(1:n*samp/2));
	#	plot(x(1:n*samp/2).*samp/n,datafrpl,"-");
	#	print -dsvg fig4.svg;

		disp(temp_inc);
		tempcor = innerpr(sigfr,datafr,0,samp/2);
		disp(tempcor), disp("");	

	endwhile
	
	tempcor=1.00;
	while(tempcor>0.97)
		temp_inc = temp_inc + inc;
		sig_new(1:n*samp +1) = x.*sin(x.*((fr_1+ temp_inc).*2*pi));
	
	#	plot(x,sig_new,"-",x,x.*0,"-");
	#	print -dsvg fig3.svg;
	
		signew_nor = sig_new./sqrt(sum(sig_new.^2));
		datafr = fft(signew_nor);
		datafrpl =abs(datafr(1:n*samp/2));
	#	plot(x(1:n*samp/2).*samp/n,datafrpl,"-");
	#	print -dsvg fig4.svg;

		disp(temp_inc);
		tempcor = innerpr(sigfr,datafr,0,samp/2);
		disp(tempcor), disp("");	

	endwhile
	


endfunction
