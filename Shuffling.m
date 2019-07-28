function Shuffling(cellsn,pairs)

%shuffling gridscore
     nshuf=10;pval=0.1;   
                                            nb=50;movm=25;asig=2;
     pgb = []; pgm = [];
     for i = 1:10% len(cellsn); 
         i
         c=cellsn(i); tic;
         pgb(i,:) = shuffleGridscoreNan(c.before,nshuf,nb,movm,asig,pval,'before');
         pgm(i,:) = shuffleGridscoreNan(c.midall,nshuf,nb,movm,asig,pval,'midall');
        toc         
     end
     
    tt = {pgb,pgm};ncls=size(pgb,1);
    cgbma = [];%zeros(len(cellsn),len(tt));
    pgbma = [];%zeros(len(cellsn),len(tt));
    for ii = 1: len(tt)
        cgbma(:,ii) = tt{ii}(:,1);
        for i = 1:ncls
            t = (tt{ii}(i,:));
            %[~,I]=sort(abs(t),'descend');if ii>6;[~,I]=sort(t,'descend');end
            [~,I]=sort(t,'descend');
            pp = find(I==1); %index of non shuffled
            pgbma(i,ii) = pp;
            if t(1)==0 %THIS IS FOR GRID SCORE CHECK
                pgbma(i,ii) = len(t);
            end
        end
    end
    %VISUALIZE
    show = [1,2];
    pp      = pgbma      <= nshuf*pval; %befire after
    pp(:,2) = pgbma(:,2) >  nshuf*(1-pval);  
    sum(pp(:,show))/ncls %LARGEST VAL (FIRST DESC)
    ppgbma = pp;
    clear pp; clear n; clear show; clear i; clear ii; clear I; clear t; clear tt;

    %pp = n*0.99 <= ptsbma; sum(pp(:,show))/len(pairs) %SMALLEST VAL (LAST DESC)
    %pp = ptsbma <= n*0.01 | n*0.99 <= ptsbma;    



%shuffling Spatial + Time
    
    par=[]; par.n=1000; par.movmean=25; par.nb = 50; t=zeros(len(pairs),par.n);
    pscb = t; pscm = t; psca = t;
    
    ptcb = t; ptcm = t; ptca = t;
    
    for i = 1:length(pairs)       
        i
        c1 = cellsn(pairs(i,1));%cells{pairs(i,1)};
        c2 = cellsn(pairs(i,2));
        %space
        tic
        pscb(i,:) = shuffleSpace2Correlations(c1.before, c2.before,par); %[pairs(i,:), ]
        pscm(i,:) = shuffleSpace2Correlations(c1.midall, c2.midall,par);
        if(length(c1.after.pt)>1 && length(c2.after.pt)>1)
            psca(i,:) = shuffleSpace2Correlations(c1.after, c2.after,par);
        else;psca(i,:) = zeros(1,par.n);end
        toc
        %time
        tic
        ptcb(i,:) = shuffleTimeCorrelations (c1.before, c2.before,par);
        ptcm(i,:) = shuffleTimeCorrelations (c1.midall, c2.midall,par);
        if(length(c1.after.pt)>1 && length(c2.after.pt)>1)
            ptca(i,:) = shuffleTimeCorrelations (c1.after, c2.after,par);
        else;ptca(i,:) = zeros(1,par.n);end 
        toc
    end
    
    return %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tt = {ptcb,ptcm,ptca,pscb,pscm,psca}%,pgcb,pgcm,pgca};
    ctsbma = zeros(len(pairs),len(tt));
    ptsbma = zeros(len(pairs),len(tt));
    for ii = 1: len(tt)
        ctsbma(:,ii) = tt{ii}(:,1);
        for i = 1:len(pairs)
            t = (tt{ii}(i,:));
            %[~,I]=sort(abs(t),'descend');if ii>6;[~,I]=sort(t,'descend');end
            [~,I]=sort(t,'descend');
            pp = find(I==1); %index of non shuffled
            ptsbma(i,ii) = pp;
            if t(1)==0 || t(1)==-2 %THIS IS FOR GRID SCORE CHECK
                ptsbma(i,ii) = len(t);
            end
        end
    end
    clear i; clear ii; clear I; clear nb; clear par; clear t; clear tt; clear pp;
    %clear ptca; clear ptcb; clear ptcm; clear psca; clear pscb; clear pscm; 
    %VISUALIZE
    n = len(ptcb); show = [1,2,3,4,5,6]; %n = len(ptcb) % ==1000
    pp = ptsbma<= n*0.01;  sum(pp(:,show))/len(pairs) %LARGEST VAL (FIRST DESC)
    pp = n*0.99 <= ptsbma; sum(pp(:,show))/len(pairs) %SMALLEST VAL (LAST DESC)
    pp = ptsbma <= n*0.01 | n*0.99 <= ptsbma;
    %pp(:,7:9) = ptsbma(:,7:9) <= n*0.01; %for GRID    
    sum(pp(:,show))/length(pairs)  %COMPARE TO SHUFF
    pptsbma = pp;
    clear pp; clear n; clear show; 
    
    
    
    
    unique(q(:,2))
    figure;histogram(q(:,2),100)
    
    
    q = [];
    for i = 1: size(pscb,1);
        t = abs(ptcb(i,:));
        [~, mi] = max(t(2:end));
        q(i,:) = [t(1), mi];
        %[s,I] = sort(abs(t),'descend');
        %pb = round(find(I==1)/length(I),2); %index of non shuffled
        %[pb pm]; pbs(end+1) = pb; pms(end+1) = pm;
    end
    unique(q(:,2))
    figure;histogram(q(:,2),100)
    
    
    %PLOT MEASUREMENTS OF GRID bef dur after
    b = []; m = []; a = [];
    for i = 1:length(pairs)
        i
        c1 = cells{pairs(i,1)};
        c2 = cells{pairs(i,2)};
        b(i,:) = apear(c1.before, c2.before); %[pairs(i,:), ]
        m(i,:) = apear(c1.midall, c2.midall);
        if(length(c1.after.pt)>1 && length(c2.after.pt)>1)
            a(i,:) = apear(c1.after, c2.after);
        end
    end
    
    h = figure(934); set(h,'Position',[-1000 10 1000 1000]);
    colstr = {'time cc smoothed at 30ms';
          'spatial cc at 0,0, no smoothing';
          'maximum gridscore';
          'gridscore of spatial cc';
         };
    stdd=3;
    t = b;              k=[1:len(t)]; k = t(:,1) < stdd*std(t(:,1));  %leaving out corr outliers
    for i = 1:3
        subplot(3,3,i)
        k = t(:,i+1)>-2 & k;
        scatter(t(k,1),t(k,i+1),'.');
        xlabel(colstr(1));ylabel(colstr(i+1));
        title('before');
    end
    t = m;              k=[1:len(t)]; k = t(:,1) < stdd*std(t(:,1));
    for i = 1:3
        subplot(3,3,3+i)
        k = t(:,i+1)>-2 & k;
        scatter(t(k,1),t(k,i+1),'.');
        xlabel(colstr(1));ylabel(colstr(i+1));
        title('during');
    end
    t = a(a(:,1)~=0,:); k=[1:len(t)]; k = t(:,1) < stdd*std(t(:,1));
    for i = 1:3
        subplot(3,3,6+i)
        k = t(:,i+1)>-2 & k;
        scatter(t(k,1),t(k,i+1),'.');
        xlabel(colstr(1));ylabel(colstr(i+1));
        title('after');
    end
    
    ax = findobj(h,'type','axes');
    for i = 1:length(ax)
        %axis(ax(i),'square');
    end
    
    
end
    
    
    
    


