%misc_paper
a=[1 1]; b=[0 1];
c=zeros(1,1000); [r,p]=ccof([c a],[c b],10)
c=zeros(1,10); [r,p]=ccof([c a],[c b],10)

%Animals in Study
t= [cellsn(cels).id]; len(unique(t))

%threshold test
bt=0.9:-0.1:0.2;
mt=0.5:-0.1:0; 
tbm=[];sbm=[];tbmp=[];sbmp=[];
for b=bt
    v=[];w=[];vp=[];wp=[]; tic
    for m=mt
        [b m]
        prs=pairsMake(cellsn,b,m);
        t=tsbma(cellsn,prs);
        [r,p]=ccof(t(:,1),t(:,2));
        v(end+1)=r;vp(end+1)=p;
        [r,p]=ccof(t(:,4),t(:,5));
        w(end+1)=r;wp(end+1)=p;
    end
    tbm(end+1,:)=v;tbmp(end+1,:)=vp;
    sbm(end+1,:)=w;sbmp(end+1,:)=wp;
    toc
end
x=repmat(bt',1,len(mt));
y=repmat(mt,len(bt),1);
figure(3); clf; set(gcf,'position',[100,100, 900,900]); 
p=(4-1)*len(bt)+5; ['thresh paper ' n2(x(p)) ' ' n2(y(p))]  % 0.5 0.2 thresh in paper  4
subplot(221);plot3(x,   y,   tbm,   'ko','markerfacecolor','g');
hold on;     plot3(x(p),y(p),tbm(p),'ko','markerfacecolor','b');
title('temporal corr by grid thresh');grid on; zlim([0 1]); 
xlabel('min grid score pre');ylabel('max grid score dur');zlabel('corr pre x dur'); 
subplot(222);plot3(x,   y,   tbmp,   'ko','markerfacecolor','g');
hold on;     plot3(x(p),y(p),tbmp(p),'ko','markerfacecolor','b');
grid on; title('temporal corr pval by grid thresh'); zlim([0 1]);
xlabel('min grid score pre');ylabel('max grid score dur');zlabel('pval pre x dur');
subplot(223);plot3(x,   y,   sbm,   'ko','markerfacecolor','m'); 
hold on;     plot3(x(p),y(p),sbm(p),'ko','markerfacecolor','b');
grid on; title('spatial corr by grid thresh'); zlim([0 1]);
xlabel('min grid score pre');ylabel('max grid score dur');zlabel('corr pre x dur'); 
subplot(224);plot3(x,   y,   sbmp,   'ko','markerfacecolor','m');
hold on;     plot3(x(p),y(p),sbmp(p),'ko','markerfacecolor','b');
grid on; title('spatial corr pval by grid thresh'); zlim([0 1]);
xlabel('min grid score pre');ylabel('max grid score dur');zlabel('pval pre x dur');

figure(5); surf(x,y,tbm); colormap jet;

%mean gridscore
t=[cellsn(cels).before]; [mean([t.gridscore]),std([t.gridscore])]
t=[cellsn(cels).midall]; [mean([t.gridscore]),std([t.gridscore])]
%percent removed
[max(percentremoved) mean(percentremoved) std(percentremoved)]*100
%recording time
t=cell2mat(cellfun(@(x) [cellsn(x(1)).midall.pt(1) cellsn(x(1)).midall.pt(end)]/60, group,'uni',0));
mean(t)
%room size
load('./data/roomdose.mat');
t=cellfun(@(x) cellfun(@(y) max(y)-min(y),x), cellfun(@(x) {x.px}, ...
    arrayfun(@(x) [cellsn(room==x).midall],unique(room),'uni',0),'uni',0),'uni',0);
cellfun(@(x) max(x),t)
figure(4), histogram(t{1},20)
figure(5), histogram(t{2},20)
figure(6), histogram(t{3},20)
figure(7), histogram(room(cels));
figure(8), histogram(room(chd))
figure(9), histogram(room(chdn));

%mean decrease in firing rate
for i = 1:len(cels)
    c = cellsn(cels(i));
    mr(i,1) = len(c.before.st)/(c.before.st(end)-c.before.st(1));
    mr(i,2) = len(c.midall.st)/(c.midall.st(end)-c.midall.st(1));
    mr(i,3) = len(c.after.st)/(c.after.st(end)-c.after.st(1));
end
100*(mean(mr(:,2))/mean(mr(:,1)))
100*(mean(mr(mr(:,3)~=inf,3))/mean(mr(:,1)))
stats(mr(:,1),1,2)
stats(mr(:,2),1,2)
stats(mr(:,3),1,2)
clear mr;

% agg shuff
par.n=1000;par.pval=0.02;par.spval=0.01;
tt = {ptcb,ptcm,ptca,pscb,pscm,psca}%,pgcb,pgcm,pgca};
    niter=size(tt{1},1);
    ctsbma = zeros(niter,len(tt));
    ptsbma = zeros(niter,len(tt));
    for ii = 1: len(tt)
        ctsbma(:,ii) = tt{ii}(:,1);
        for i = 1:niter
            t = (tt{ii}(i,:));
            [~,I]=sort(abs(t),'descend'); %if ii>6;[~,I]=sort(t,'descend');end
            %[~,I]=sort(t,'descend');
            pp = find(I==1); %index of non shuffled
            ptsbma(i,ii) = pp;
            if t(1)==0 || t(1)==-2 %THIS IS FOR GRID SCORE CHECK
                %ptsbma(i,ii) = len(t);
            end
        end
    end
clear i; clear ii; clear I; clear t; clear tt; clear pp;
'(+) corr'
100*sum(ptsbma<= par.n*par.spval&ctsbma>0)/len(pairs)
'(-) corr'
100*sum(ptsbma<= par.n*par.spval&ctsbma<0)/len(pairs)
%pptsbma = pp;
pptsbma= ptsbma <= par.n*par.spval;
'corr'
100*sum(pptsbma)/len(pairs)
[~,~,s] = ccof(ctsbma(:,1),ctsbma(:,2),5)
[~,~,s] = ccof(ctsbma(:,1),ctsbma(:,3),5)
[~,~,s] = ccof(ctsbma(:,4),ctsbma(:,5),5)
[~,~,s] = ccof(ctsbma(:,4),ctsbma(:,6),5)

%out group corrs
load('.data/outgroupcorrs_v2');
p.nb = 50; p.movmean=25; csb=[];csm=[];ctb=[];ctm=[];cii=[];
for i = 1:len(cels)
    fprintf(' %d',i)
    tic
    c1 = cellsn(cels(i));
    for j = i+1:len(cels)
        c2 = cellsn(cels(j));
        if ~isequal(sprintf('%s%s',c1.date,c1.id),sprintf('%s%s',c2.date,c2.id));
            csb(end+1) = spaceCorrelationSmoothed(c1.before,c2.before,p);
            csm(end+1) = spaceCorrelationSmoothed(c1.midall,c2.midall,p);
             [t1, t2] = createMsSpikeTrain(c1.before.st, 100, c2.before.st);
             ctb(end + 1) = timeCorrelationSmoothed(t1,t2,p);
             [t1, t2] = createMsSpikeTrain(c1.midall.st, 100, c2.midall.st);
             ctm(end + 1) = timeCorrelationSmoothed(t1,t2,p);
            cii(end+1,:) = [c1.ind,c2.ind];
        end
    end
    toc;
end

%[mean(abs(ctb)), mean(abs(ctm))] %OUT OF GROUP
%[mean(abs(csb)), mean(abs(csm))] %OUT OF GROUP
[~,~,s]=ccof(ctb,ctm,5)
[~,~,s]=ccof(csb,csm,5)
clear cii; clear c1; clear c2; clear i; clear j;

%time vs space corr
[~,~,s] = ccof(ctsbma(:,1),ctsbma(:,4),5)
[~,~,s] = ccof(ctsbma(:,2),ctsbma(:,5),5)
[~,~,s] = ccof(ctsbma(:,3),ctsbma(:,6),5)


% random hd clusters
trandclusters = [];
n=100000;tic %13sec for 1e5 | 170sec for 1e6 | 107 pairs
for i = 1:n
    rphd =  logical(zeros(len(pairs),1));rphdn = rphd;
    t=randperm(len(pairs),len(phd) + len(phdn))'; 
    rphd(t(1:len(phd)))=true;rphdn(t(len(phd)+1:end))=true;    
    trandclusters(i,:)=[ccof(ctsbma(rphd,1),ctsbma(rphd,2)),ccof(ctsbma(rphdn,1),ctsbma(rphdn,2)),...
                        ccof(ctsbma(rphd,4),ctsbma(rphd,5)),ccof(ctsbma(rphdn,4),ctsbma(rphdn,5))];
end;toc
t12 =  abs(ccof(ctsbma(phd,1),ctsbma(phd,2)) - ccof(ctsbma(phdn,1),ctsbma(phdn,2)));
tt12 = abs(trandclusters(:,1)-trandclusters(:,2));
t34 =  abs(ccof(ctsbma(phd,4),ctsbma(phd,5)) - ccof(ctsbma(phdn,4),ctsbma(phdn,5)));
tt34 = abs(trandclusters(:,3)-trandclusters(:,4));
100*sum(abs(t12)>abs(tt12))/len(trandclusters)
100*sum(abs(t34)>abs(tt34))/len(trandclusters)
%stats(tt12);stats(tt34);
figure(12098);
subplot(211);histogram(tt12,100);title('time');
line([t12, t12], ylim, 'LineWidth', 2, 'Color', 'r');
subplot(212);histogram(tt34,100);title('spatial');
line([t34, t34], ylim, 'LineWidth', 2, 'Color', 'r');
suptitle(sprintf('histogram of n=%d runs of random clusters \n[nonHD - HD] of corr(before,during)',10000)) 
clear t12;clear t34; clear tt12;clear tt34; clear rphd; clear rphdn;



% Average Angle
rst=0.4; t = [cellsn.before]; tl=arrayfun(@(z) len(z.st),t); ta = [t.('rayleigh_angle')]; 
ts = [t.('rayleigh_score')]; tis = ts>rst&tl>=100; ta = ta(tis); ta=toCol(ta);
rad2deg([circ_mean(ta) circ_std(ta)])
rst=0.4; t = [cellsn.midall]; tl=arrayfun(@(z) len(z.st),t); ta = [t.('rayleigh_angle')]; 
ts = [t.('rayleigh_score')]; tis = ts>rst&tl>=100; ta = ta(tis); ta=toCol(ta);
rad2deg([circ_mean(ta) circ_std(ta)])
rst=0.4; t = [cellsn(chd).midall]; tl=arrayfun(@(z) len(z.st),t); ta = [t.('rayleigh_angle')]; 
ts = [t.('rayleigh_score')]; tis = ts>rst&tl>=100; ta = ta(tis); ta=toCol(ta);
rad2deg([circ_mean(ta) circ_std(ta)])


