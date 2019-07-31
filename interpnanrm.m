%interpolates nans in ratemap
function rmintrp=intrepnanrm(rm)
    [m n]=size(rm);
    [x y]=meshgrid(1:n,1:m);
    x=x(:);y=y(:);rm=rm(:);
    nni=isnan(rm); nx=x;ny=y;
    nx(nni)=[];ny(nni)=[];rm(nni)=[];
    f = scatteredInterpolant(nx,ny,rm);
    rmintrp=f(x,y);
    rmintrp=reshape(rmintrp,m,n);
end


%{
rm=cellsn(2).before.rm;
%rm=rm(25:50,:)
figure(1);
subplot(211);
rmt=rm; %rmt(isnan(rmt))=99;
imgsc(rmt)
subplot(212);
imgsc(interpnanrm(rm),3);
%}
