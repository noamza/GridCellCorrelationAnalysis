function f = plotfft(ax,x,fs)
    y=fft(x-mean(x));
    l = length(x);
    f = fs*(1:floor(l/2))/l; %/2??
    y(1:1000)=0; %deletes first few;
    p = plot(ax,f,abs(y(1:floor(l/2))));
    %plot(ax,f,[toCol(real(y(1:floor(l/2)))) , toCol(complex( y(1:floor(l/2) ) ) ) ] );
    hold(ax,'on'); plot(ax,[7.9999 8],[ylim(ax)],'r--','linewidth',1.2); hold(ax,'off'); %max(y(1:round(l/2)) )
    xlabel(ax,'Hz');
end