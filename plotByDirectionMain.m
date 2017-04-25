function [h, vers] =  plotByDirectionMain(params)
    %PARAMETERS
    params.gid
    binl = params.binspike;%0.06 0.1;%.3
    lag = params.lag;
    nbins = round(360/params.bindeg);
    fprintf('nbins %d\n', nbins);
    sigma = params.sigma;
    v = params.version;
    params.number_degree_bins = nbins; params.lag_max_secs = lag; params.time_bin_secs = binl;
    k = fieldnames(params.groups); %dt = 0.02;
    g =  params.groups.(k{params.gid});
    params.grid_thresh;%0.5; % was .7
    sesh = params.sesh;
    parent = params.parent;
    %fprintf('%s %d %d\n', k{l}, l, length(groups.(k{l})));
    
    if isfield(params,'fig')
        h = params.fig;
    else
        h = figure('Position', [0, 0, 2000, 1000]);%(2*m+5)*120 , 140*r + height]); %
        set(gca,'LooseInset', get(gca,'TightInset'));
    end;
    good = [g(1)]; bad = [g(1)]; %find low gridscore cells
    for j = 1:length(g) %put low grid scores at end
        if g(j).before.gridscore > params.grid_thresh;
            good(end+1) = g(j); %#ok<*AGROW>
        else bad(end+1) = g(j);
        end; 
    end
    bad = bad(2:end); good = good(2:end); %[good(2:end) bad(2:end)]; Only show good cells
    if isfield(params,'g')
%         good =  params.groups.(k{params.gid});
    else
    %%% plot each cell %%%
    c = length(good); r = c;  if r == 1; r=3; end;
    if ~isempty(good) && length(good) >= 2        
        
        colormap jet;
        %for each cell, for each direction x each other cell
        for i = 1:length(good)-1
            for j = i+1:length(good)
                
                if strcmp(params.sesh,'before')
                    c1 = good(i).before; c2 = good(j).before; clast = good(r).before;
                elseif strcmp(params.sesh,'after')
                    c1 = good(i).after; c2 = good(j).after;   clast = good(r).after;
                else
                    %ADD CASE FOR BAD MIDDLE
                    c1 = good(i).middle{params.sesh}; c2 = good(j).middle{params.sesh};
                    clast = good(r).middle{params.sesh};
                    sesh = sprintf('mid%d',params.sesh);
                end 
                c1.ind = good(i).ind; c2.ind = good(j).ind;  clast.ind = good(r).ind;                
                [p co] = compareByMovingDirection(c1, ...
                    c2, params); %j
                subplot(r, c, c*(i-1) + j-1,'Parent', parent);
                %X = [0:length(p{1})-1]'-(length(p{1})-1)/2; X=X*params.time_bin_secs;
                X = p{1,2}';
                Y = reshape(cell2mat(p(:,1)),[],nbins);
                %[~,ax,AX] = plotmatrix(X,Y,'-');
                Y2 = Y'; ax = gca; %AX = gca;
                %Plot correlation by degree figs
                %imagesc(imgaussfilt(Y2, 2));set(gca,'ydir','normal'); colormap jet;  
                %any(isnan(Y2(:)) useful
                Y3 = imgaussfilt([Y2;Y2;Y2], sigma,'FilterDomain','spatial'); %for nans
                if any(~isnan(Y3(:))) 
                    Y3 = Y3(size(Y2,1)+1:size(Y2,1)*2,1:size(Y2,2));
                    imagesc(Y3);set(gca,'ydir','normal'); 
                    cmax = max(abs(Y3(:)));
                else %NAN problem
                    imagesc(Y2);
                    cmax = max(abs(Y2(:)));
                end
                colormap jet; set(gca,'ydir','normal');caxis([-cmax cmax]);
                
                for q=1:length(ax)
                    %axis(ax(q),'off');
                    axis(ax(q),'tight')
                    %ylim(ax(q),[min(Y(:)) max(Y(:))]);
                    %ylabel(ax(q),sprintf('%.0f°',mean(p{q,3}))); %mean of bin edges (midway)
                    set(ax(q),'xtick',[1 round(length(X)/2) length(X)]);
                    set(ax(q),'xticklabel',[X(1) 0 X(end)]);%for imagesc
                    %set(ax(q),'ytick',[]);
                    %set(ax(q),'yticklabel',[]);
                    set(ax(q),'ytick',1:length(p(:,3)));
                    set(ax(q),'yticklabel',mean(cell2mat(p(:,3))')); %#ok<*UDIM>
                    set(ax(q),'FontSize',6);%axis(ax(q),'equal');%axis(ax(q),'off');%title(ax(i),{'Very Nice'})
                end
                %}
                if i==1 && j==2 %used to be AX
                    axis(ax,'tight')
                    %set(AX,'FontSize',8)
                    set(ax,'fontweight','bold');
                    ylabel(ax,sprintf('MD bin (b=%.1f°)',360/nbins));
                    xlabel(ax,sprintf('%s sbin%.2fs','Xcorr lag(s)', ...
                        params.time_bin_secs));
                end
                title(sprintf('c%d*c%d |mxcrr|%.2f',c1.ind,c2.ind,round(max(abs(Y(:))),2)));
                %CROSS RMAT
                subplot(r, c, r*r-(c*(i-1) + j-1),'Parent', parent);
                cc = xcorr2(c1.rm,c2.rm);
                imagesc(cc);xlabel(sprintf('c%d x c%d',c1.ind,c2.ind)); %upside down?
                axis equal; axis tight; set(gca,'ydir','normal','xticklabel',[],'yticklabel',[]);colormap jet;
                hold on; plot(size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7)
                
            end
            %plot last rate mat ?
            subplot(r,c,i*r,'Parent', parent); imagesc(c1.rm);
            title(sprintf('c%d(gsc%.1f)',c1.ind,round(c1.gridscore,1)));
            axis equal; axis tight; set(gca,'ydir','normal','xticklabel',[],'yticklabel',[]);colormap jet;
        end
        %plot last rate mat ?
        subplot(r,c,r*r,'Parent', parent); imagesc(clast.rm);
        title(sprintf('c%d(gsc%.1f)',clast.ind,round(clast.gridscore,1)));
        axis equal; axis tight; set(gca,'ydir','normal','xticklabel',[],'yticklabel',[]);colormap jet;        
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

