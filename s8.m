%based on 'Revisopm.m'
function s8()

    a=[],b=[],c=[],d=[],e=[],f=[],m=[],n=[],x=[],y=[],z=[];
    tic
    %load('C:\Noam\Data\muscimol\aclls15min','aclls'); %personal
    load('.\data\pairs','pairs');
%     load('.\\data\\dxdyrate','dxywinrd','dxywinrn','gridwin','autowin');
    load('.\\data\\dxdyrate15','dxywinrd15');
    %cels = unique([pairs(:,1),pairs(:,2)])';
    toc;

    nbins = 100; twins = [1 2 3 5 10] ; %window in secs
    s = {'before', 'midall','after'}; ts = {'pre','dur','post'};

    fig = figure(81); clf;  fs = 12;
    set(fig,'Color','w', 'Position', [600 0 1200 800]);
    g0 = uix.GridFlex('Parent',fig,'Spacing',5, 'BackgroundColor','w');
    gtop = uix.GridFlex('Parent',g0,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
    axes('Parent',uicontainer('Parent',gtop,'BackgroundColor','w'),'visible','off');
    text(0,0.5,'A','fontweight','bold','fontsize',fs);
    gA = uix.GridFlex('Parent',gtop,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');    
    gbot = uix.GridFlex('Parent',g0,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');
    axes('Parent',uicontainer('Parent',gbot,'BackgroundColor','w'),'visible','off');
    text(0,0.5,'B','fontweight','bold','fontsize',fs);
    gB = uix.GridFlex('Parent',gbot,'Spacing',5, 'BackgroundColor','w','DividerMarkings','off');    
    %A - imgsc 
    a = [12 13 16]; d = 16; d = 50-d:51+d; %b=1;
    axg = [];  twin=1; dxyr=dxywinrd15;
    for b = a
        %c = aclls(pairs(b,1)); d = aclls(pairs(b,2));
        x=1; c = dxyr{b,twin}.(s{x}); c=c(d,d);
        axg(end+1) = axes('Parent',uicontainer('Parent',gA));
        axg(end)= imgsc(c,2); title(sprintf('p%d 1s %s',b,ts{x})); axis off;
        x=2; c = dxyr{b,twin}.(s{x}); c=c(d,d);
        axg(end+1) = axes('Parent',uicontainer('Parent',gA));
        axg(end)= imgsc(c,2); title(sprintf('p%d 1s %s',b,ts{x}),'color','red'); axis off;
    end
    set(gA,'Widths', zeros(1,len(a)*2)-1, 'Heights', -1)
    set(gtop,'Widths', -1, 'Heights', [ 25, -1]);
    %end gtop A
    % gbottom SCATTER
    bxg=[]; m=1;w=2; mm=['corr ' ts{m} ' ' ts{w}];
    for h =1:len(twins)
        twin=twins(h);
        g={}; %scatter
        for w = 2:3 %scatter
            a=[];
            for i = 1:len(dxyr)
                %tic; [twin i]
                %corr sesh vs sesh for dxdy
                t=corr(dxyr{i,twin}.(s{m})(:),dxyr{i,twin}.(s{w})(:),'rows','complete');
                a(end+1,:)=t; %toc
            end
            %g(:,h)=a; %g{h}=a(~isnan(a));%
            g{w}=a; %scatter
            %gridwin{twin}=a;
        end
        bxg(end+1) = axes('Parent',uicontainer('Parent',gB));
        y=toCol(g{2});x=toCol(g{3});t=~isnan(y)&~isnan(x);x=x(t);y=y(t); bxg(end)=scatter(x,y); hold on;
        f = fit(x, y,'poly1');plot(x,f(x),'-'); axis('tight'); axis square; [a b]=ccof(x,y);
        text(0.1,0.9,sprintf('a=%.2f r=%.2f p=%.3f',round(f.p1,2), round(a,2), round(b,3)),'Units','normalized');
        title([n2(twin) 's']); xlabel('pre vs post corr'); ylabel('pre vs dur corr'); box on
        e=slim(gca);xlim(e);ylim(e);
    end
    %scatter
    t=['correlation of dxdy by time window ' z];
    %suptitle(t);
    set(gB,'Widths', [-1 -1 -1], 'Heights', [-1 -1])
    set(gbot,'Widths', -1, 'Heights', [ 25, -1]);
    set(g0,'Widths', -1, 'Heights', [ -1, -4]);
    a = findobj(fig,'type','UIContainer');
    for i = 1:len(a)
        a(i).BackgroundColor = 'w';
    end
    %saveas(f,['./figs/scatter ' t '.png']);
   
    
            % % % % graphics imagesc 
        d = 16; d = 50-d:51+d; %b=1;
        pr = [200 10 1800 900];  t=2; e = 2;% f=5;f=ceil(sqrt(f));%f-#pairs %t 
        for twin = 1%twins
            figure(1000+10*twin+t); clf; set(gcf,'position',pr);
            for i = 1:20%(pairs)
                subplot(4,e*5,e*i-e+1);x = 1; %PRE 1 DUR 2
                a = dxywinrd15{i,twin}.(s{x}); a=a(d,d); %b = gridscore2(a,t);
                imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                %imagesc(xcorr2(imgaussfilt(a, t,'FilterDomain','spatial')));
                %a(isnan(a))=-1;%imagesc(a);
                colormap jet; axis off; axis square;title(sprintf('p%d:%s',i,ts{x}));
                
                subplot(4,e*5,e*i); x = 2;%PRE 1 DUR 2
                a = dxywinrd15{i,twin}.(s{x}); a=a(d,d); % b = gridscore2((a),t);
                imagesc(imgaussfilt(a, t,'FilterDomain','spatial'));
                %imagesc(imgaussfilt(xcorr2(a), t,'FilterDomain','spatial'));
                colormap jet; axis off; axis square;title(sprintf('p%d:%s',i,ts{x}),'color','red');
            end
            %suptitle(sprintf('%s %ds: dxdy-spike/dxdy-pos and gridscore by pair sigma=%d',s{sr},twin,t));
            suptitle(sprintf('%ds: dxdy-spike/dxdy-pos and gridscore by pair sigma=%d',twin,t));
            %saveas(gcf, sprintf('./figs/drift_diff_%ds_sig%d_zoom_dmet.png',twin,t) );
        end
    
    
    %{
    t=[aclls.(s{1})];d=cellfun(@(x) max(diff(x)),{t.pt},'uni',1);
    t=[aclls.(s{2})];e=cellfun(@(x) max(diff(x)),{t.pt},'uni',1); off=max([d e]);%0.04
    dxywinrd15={};
    for h = 1:len(twins)
        twin = twins(h);
        for j=1:len(pairs)
            ['a ' n2(twin) ' ' n2(j)]
            tic
            for ss= {'midall'}%[s(:)']
                a = aclls(pairs(j,1)).(ss{:}); b = aclls(pairs(j,2)).(ss{:});
                rmsa = nan(nbins); rmta=rmsa; tssi = 1;tppi = 1; %norm=zeros(nbins);
                for i = 1:len(a.st)
                    sr = driftwin(twin,a.st(i),a.sx(i),a.sy(i),b.st,b.sx,b.sy,tssi,0);
                    pr = driftwin(twin,a.st(i),a.sx(i),a.sy(i),b.pt,b.px,b.py,tppi,off);
                    %z=~ismember(sr.a(:),pr.a(:)); assert(sum(z)==0);
                    rms = driftmap(sr.a,max(a.px),max(a.py)); tssi = sr.i;
                    rmt = driftmap(pr.a,max(a.px),max(a.py)); tppi = pr.i;
                    %there should not be bins with spike and no time.
                    rmt(~isnan(rms) & isnan(rmt)) =rms(~isnan(rms) & isnan(rmt)); %ADD BACK
                    %should not be bins with more spikes than time?
                    rmt(rms>rmt)=rms(rms>rmt);
                    t = nansum(cat(nbins,rmsa,rms),nbins);
                    t(isnan(rmsa) & isnan(rms))=nan; %only maintain nan if both values are not nan.
                    rmsa=t;
                    t = nansum(cat(nbins,rmta,rmt),nbins); %sets all nans to 0's
                    t(isnan(rmta) & isnan(rmt))=nan; %only keep nan if both values are nan.
                    rmta=t;
                end

                %assert(isequal(norm==0,isnan(dxyr)));
                dxyr=rmsa./rmta;
                %assert( max(dxyr(:) ) <= 1 );
                dxywinrd15{j,twin}.(ss{:}) = dxyr;
            end
            toc;
        end
    end
    for h = 1:len(twins)
        twin = twins(h);
        for j=1:len(pairs)
            ss= 'after';
                dxywinrd15{j,twin}.(ss) = dxywinrd{j,twin}.(ss);
        end
    end
    
    %drift window
    %has rounding error between discrete x y vals and time x y
    function r = driftwin(twin,t,tx,ty,st,sx,sy,tsi,off)
    twin=twin+off;
    dxya = []; tri = tsi;
    emd=length(st); tmin=t-twin; tmax=t+twin;
    while tsi <= emd
        if   st(tsi) < tmin  % st(tsi) + twin < t
            tri=tsi;%?????
        elseif  tmax < st(tsi) %t < st(tsi)-twin
            tsi = emd;
        else
            %if b.st(k)-twin <= t && t <= b.st(k)+twin
            dx =  sx(tsi) - tx; %no need to round.
            dy =  sy(tsi) - ty; %no need to round.
            dxya(end+1,:) = [dx dy];
        end
        tsi=tsi+1;
    end
    r.a=dxya; r.i = tri;
    %assert (tri~=-1)
    end

    %rate map
    function rm = driftmap(dxya,mpx,mpy) %a contains all diffs of x and y;
    rm = zeros(nbins);
    if length(dxya)>2
        dx = dxya(:,1);dy = dxya(:,2); %px = px+max(px); py = py+max(py);
        binx=-mpx:2*mpx/nbins:mpx;
        biny=-mpy:2*mpy/nbins:mpy;
        pxi = discretize(dx, binx);
        pyi = discretize(dy, biny);
        rm = rm + accumarray([pyi, pxi], 1, [nbins nbins]);
    else
        'empty drift map';
    end
    %rm(rm==0)=nan;
    end
    %}

end
