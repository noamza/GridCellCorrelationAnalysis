%misc_paper


%mean decrease in firing rate
for i = 1:len(cels)
    c = cellsn(cels(i));
    mr(i,1) = len(c.before.st)/(c.before.st(end)-c.before.st(1));
    mr(i,2) = len(c.midall.st)/(c.midall.st(end)-c.midall.st(1));
    mr(i,3) = len(c.after.st)/(c.after.st(end)-c.after.st(1));
end
100*(mean(mr(:,2))/mean(mr(:,1)))
100*(mean(mr(mr(:,3)~=inf,3))/mean(mr(:,1)))
clear mr;

% agg shuff
par.n=1000;par.pval=0.02;par.spval=0.01;  
'(+) corr'
sum(ptsbma<= par.n*par.spval&ctsbma>0)/len(pairs)
'(-) corr'
sum(ptsbma<= par.n*par.spval&ctsbma<0)/len(pairs)
%pptsbma = pp;
pptsbma= ptsbma <= par.n*par.spval;
'corr'
sum(pptsbma)/len(pairs)

[~,~,s] = ccof(ctsbma(:,1),ctsbma(:,2),5)
[~,~,s] = ccof(ctsbma(:,1),ctsbma(:,3),5)

[~,~,s] = ccof(ctsbma(:,4),ctsbma(:,5),5)
[~,~,s] = ccof(ctsbma(:,4),ctsbma(:,6),5)

