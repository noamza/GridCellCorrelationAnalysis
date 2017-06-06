function [h, vers] =  plotByDirectionMain(params)
    %PARAMETERS
    params.gid
    binl = params.binspike;%0.06 0.1;%.3
    lag = params.lag;
    nbins = round(360/params.bindeg); %rounds good! 
    fprintf('nbins %d\n', nbins);
    sigma = params.sigma;
    v = params.sesh;
    params.number_degree_bins = nbins; params.lag_max_secs = lag; params.time_bin_secs = binl;
    %k = fieldnames(params.groups); %dt = 0.02;
    g =  params.groups{params.gid};
    params.grid_thresh;%0.5; % was .7
    sesh = params.sesh
    if isa(params,'App')
        h = params.UIFigure;
        parent = params.parent;
    elseif isfield(params,'parent')
        parent = params.parent;
        h = params.fig;
    else
        h = figure('Position', [0, 0, 2000, 1000]);%(2*m+5)*120 , 140*r + height]); %
        set(gca,'LooseInset', get(gca,'TightInset')); 
        parent = h; %CHECK <<
    end
    good = [g(1)]; bad = [g(1)]; %find low gridscore cells
    for j = 1:length(g) %put low grid scores at end
        if g(j).before.gridscore > params.grid_thresh;
            good(end+1) = g(j); 
        else
            bad(end+1) = g(j);
        end 
    end
    bad = bad(2:end); good = good(2:end); %[good(2:end) bad(2:end)]; Only show good cells
    m.good = good;
    %%% plot each cell %%%
    c = length(good); r = c;  if r == 1; r=3; end; z = 0;
    if ~isempty(good) && length(good) >= 2        
        %wt = parent.Position(3); dwt = floor(wt/c);
        %ht = parent.Position(4); dht = floor(ht/c);
        %for each cell, for each direction x each other cell
        mxcrs = [];
        corrAxes = [];
        for i = 1:length(good)-1
            for j = i+1:length(good)
                z = z + 1;
                if strcmp(params.sesh,'before')
                    c1 = good(i).before; c2 = good(j).before; clast = good(r).before;
                elseif strcmp(params.sesh,'midall')
                    c1 = good(i).midall; c2 = good(j).midall;   clast = good(r).midall;
                elseif strcmp(params.sesh,'after')
                    c1 = good(i).after; c2 = good(j).after;   clast = good(r).after;
                else
                    %ADD CASE FOR BAD MIDDLE
                    c1 = good(i).middle{params.sesh}; c2 = good(j).middle{params.sesh};
                    clast = good(r).middle{params.sesh};
                    sesh = sprintf('mid%d',params.sesh);
                end 
                c1.ind = good(i).ind; c2.ind = good(j).ind;  clast.ind = good(r).ind;                
                [p co] = compareByMovingDirection(c1,c2, params);
                
                az = subplot(r, c, c*(i-1) + j-1,'Parent', parent);
                %X = [0:length(p{1})-1]'-(length(p{1})-1)/2; X=X*params.time_bin_secs;
                X = p{1,2}';
                if nbins == 1
                    Y = p{1,1}';
                    rawmax = max(abs(Y(:)));
                    if sigma ~= 0
                       gwindow = 30; %gaussian win: 15 vs 30 very similar default matlab is 5;
                        Y = smoothts(Y,'g',gwindow,sigma);
                    end
                    plot(az, X, Y,'linewidth', 3); axis(az, 'tight'); ylim(az, [-1 1]); %corr not normalized
                    az.Color = 'black';
                    
                else
                    Y = reshape(cell2mat(p(:,1)),[],nbins);
                    %[~,ax,AX] = plotmatrix(X,Y,'-');
                    Y = Y'; %Y2 = Y' %imagesc(imgaussfilt(Y2, 2));set(gca,'ydir','normal');
                    rawmax = max(abs(Y(:)));
                    if sigma ~= 0
                        Ysm = imgaussfilt([Y;Y;Y], sigma,'FilterDomain','spatial'); %for nans
                        Ysm = Ysm(size(Y,1)+1:size(Y,1)*2,1:size(Y,2));
                        if any(~isnan(Ysm(:))) %if any are not nan
                            Y = Ysm;
                        end    
                    end
                    imagesc(az, Y); cmax = max(abs(Y(:))); colormap(az,'jet'); 
                    set(az,'ydir','normal'); caxis(az, [-cmax cmax]);
                    axis(az,'tight')
                    set(az,'xtick',[1 round(length(X)/2) length(X)]);
                    set(az,'xticklabel',[X(1) 0 X(end)]);%for imagesc
                    set(az,'ytick',1:length(p(:,3)));
                    set(az,'yticklabel', mean(cell2mat(p(:,3))')); %#ok<*UDIM>
                    if i==1 && j==2 %used to be AX)
                        set(az,'fontweight','bold');
                        ylabel(az,sprintf('MD bin (b=%.1f°)',360/nbins));
                        xlabel(az,sprintf('%s sbin%.2fs','Xcorr lag(s)', ...
                        params.time_bin_secs));
                    end
                end
                mxcrs = [mxcrs rawmax];
                corrAxes = [corrAxes az];
                title(az, sprintf('%dx%d|mx|%.2f',c1.ind,c2.ind,rawmax));
                set(az,'FontSize',7);
                %loop for plotmatrix(?) for i=1:length(ax) %title(ax(i),{'Very Nice'}) axis(ax(i),'off');        
                %}
                %CROSS RMAT
                az = subplot(r, c, r*r-(c*(i-1) + j-1),'Parent', parent); %axis(az);
                cc = xcorr2(c1.rm,c2.rm);
                imagesc(az, cc); xlabel(az, sprintf('c%d x c%d',c1.ind,c2.ind)); %upside down?
                axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
                colormap(az,'jet');
                hold(az,'on'); plot(az, size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7)
                
            end
            %plot rate mat 
            az = subplot(r,c,i*r,'Parent', parent); imagesc(az, c1.rm);
            title(az, sprintf('c%d(gsc%.1f)',c1.ind,round(c1.gridscore,1)));
            axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
            colormap(az,'jet');
        end
        %ADJUST MXCORRS
        mxcr = max(abs(mxcrs));
        if nbins == 1 
            for i = 1:length(corrAxes)
                ylim(corrAxes(i), [-mxcr mxcr]);
            end
        end
        %plot last rate mat ?
        az = subplot(r,c,r*r,'Parent', parent); imagesc(az, clast.rm);
        title(az, sprintf('c%d(gsc%.1f)',clast.ind,round(clast.gridscore,1)));
        axis(az,'equal'); axis(az,'tight'); set(az,'ydir','normal','xticklabel',[],'yticklabel',[]);
        colormap(az,'jet');        
        gid = params.gid;
        vers = sprintf('%s_G%d_%s_deg%d_lag%.1f_blen%.2f_sig%d_',...
            v,gid,sesh,round(360/nbins),lag,binl,sigma);
        titl = sprintf('%s_moving_direction_xcorr_ncells%d',vers,length(good));
%         set(suptitle(titl),'Interpreter', 'none'); %PUT SUP TITLE AFTER ALL SUBPLOT COMMANDS << test
%         handles.fig = h;
%         handles.ax = gca;
    else
        vers = '';
        close(h);
    end
end

