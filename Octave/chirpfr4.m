function chirpfr4(fr,n,samp,fr_end,t_n)
	#plot the correlation values of template and signal at various phases
	#noise added is 0.1*max of the template
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
	t_x=linspace(0,t_n, (t_n*samp)+1);
	fr_1=(t_x.*k)+fr;
	sig(1:t_n*samp +1) = t_x.*sin(t_x.*(fr_1.*2*pi));
	plot(x,sig,"-",t_x,t_x.*0,"-");	
	print -dsvg fig1.svg;

	#noise
	data = (zeros(1,n*samp +1));
	datawn = rand(1, n*samp +1);
	datawn = (datawn-0.5)*0.1*max(abs(sig));

	
	plot(x,datawn,"-",t_x,t_x.*0,"-");
	print -dsvg fig2.svg;

	datawn(1:(t_n*samp)+1) = datawn(1:(t_n*samp)+1)+sig(1:t_n*samp +1);
	plot(x,datawn);
	print -dsvg fig3.svg;

	datafr1= datawn./sqrt(sum(datawn.^2));
	sigfr1 = sig./sqrt(sum(sig.^2));
	
	sigfr = fft(sigfr1);
	sigfrpl =abs(sigfr(1:n*samp/2));
	plot(x(1:n*samp/2).*samp/n,sigfrpl,"-");
	print -dsvg fig4.svg;

	
	datafr = fft(datafr1);
	datafrpl =abs(datafr(1:n*samp/2));
	plot(x(1:n*samp/2).*samp/n,datafrpl,"-");
	print -dsvg fig5.svg;
	
	sig_sig_innrpr = innerpr(sigfr,sigfr,0,samp/2)
	sig_data_innrpr = innerpr(sigfr,datafr,0,samp/2)

	t_cor=cor(datawn,sig)
	f_cor=cor(datafrpl,sigfrpl)
#	cor(datawnfrpl,sigfrpl)
	
#	cor(sigfrpl,datafrpl)


endfunction
