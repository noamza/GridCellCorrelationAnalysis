function [m, h] =  driftPlot(m)
    t = '';
    if isempty(m.windows) %~isfield(m, 'windows')
        m.windows = {}; %delete
        for i = 1:length(m.g)
            if strcmp(m.sesh,'before')
                s = m.g(i).before;
            elseif strcmp(m.sesh,'midall')
                s = m.g(i).midall;
            elseif strcmp(m.sesh,'after')
                s = m.g(i).after;
            else
                s = m.g(i).middle{m.sesh};
                t = ['mid' num2str(m.sesh)];
            end
            %   driftWindow()
            m.windows{i} = driftWindow(s, m.winsecs);
        end
        m.crossed = cell(length(m.windows{1})+1,length(m.g)-1,length(m.g));
        for e = 1: length(m.windows{1})
            for i = 1:length(m.g)-1
                for j = i+1:length(m.g)
                    m.crossed{e,i,j} = xcorr2(m.windows{i}{e}.rm, m.windows{j}{e}.rm);
                    %m.crossed{e,i,j} = Cross_Correlation(m.windows{i}{e}.rm, m.windows{j}{e}.rm);
                    if sum(size(m.crossed{end,i,j}) == 0) %initialize
                        m.crossed{end,i,j} = m.crossed{e,i,j};
                    else
                        %[e i j; size(m.crossed{end,i,j}) 0 ; size(m.crossed{e,i,j}) 0]
                        mins = min( size(m.crossed{end,i,j}), size(m.crossed{e,i,j}) ); %amazingly this works
                        %running average ( previous mean * (count -1)') + new value ) / count
                        m.crossed{end,i,j} = (m.crossed{end,i,j}...
                            (1:mins(1),1:mins(2))*(e-1) +  m.crossed{e,i,j}(1:mins(1),1:mins(2)))/e;
                    end
                    
                end
            end
        end
    end
    %title
    if m.wini <= length(m.windows{1})
        m.title = ['Session-' num2str(m.sesh) sprintf(' [%ds - %ds]', ...
            round(m.windows{1}{m.wini}.start),...
            round(m.windows{1}{m.wini}.end))];
    else
        m.title = ['Session-' num2str(m.sesh) ' Mean Xcorr | Complete RM)'];
    end
    %plot it up
    h = plotFrame(m);
    
end

function h = plotFrame(m)
    
    if ~isfield(m,'parent')
        m.parent = figure('Position', [0, 0, 2000, 1000]);%(2*m+5)*120 , 140*r + height]); %
        set(gca,'LooseInset', get(gca,'TightInset'));
    end
    h = m.parent;
    c = length(m.g); r = c;  if r == 1; r=3; end; z = 0;
    
    %grid scatter
    a = [length(m.g)]; b = [length(m.g)];
    for i = 1:length(m.g)
        a(i) = m.g(i).before.gridscore;
        b(i) = m.g(i).midall.gridscore;
    end
    
    if ~isempty(m.g) && length(m.g) >= 2
        for i = 1:length(m.g)-1
            for j = i+1:length(m.g)
                z = z + 1;
                if strcmp(m.sesh,'before')
                    c1 = m.g(i).before; c2 = m.g(j).before; clast = m.g(r).before;
                elseif strcmp(m.sesh,'midall')
                    c1 = m.g(i).midall; c2 = m.g(j).midall; clast = m.g(r).midall;
                elseif strcmp(m.sesh,'after')
                    c1 = m.g(i).after; c2 = m.g(j).after; clast = m.g(r).after;
                else
                    %ADD CASE FOR BAD MIDDLE
                    c1 = m.g(i).middle{m.sesh}; c2 = m.g(j).middle{m.sesh};
                    clast = m.g(r).middle{m.sesh};
                end
                c1.ind = m.g(i).ind; c2.ind = m.g(j).ind; clast.ind =  m.g(r).ind;
                w1.ind = c1.ind; w2.ind = c2.ind; wlast.ind = clast.ind;
                if m.wini <= length(m.windows{1})
                    w1.rm    = m.windows{i}{m.wini}.rm; %m.wini
                    w2.rm    = m.windows{j}{m.wini}.rm;
                    wlast.rm = m.windows{r}{m.wini}.rm;
                else %plot aggragate
                    w1.rm    = c1.rm;
                    w2.rm    = c2.rm;
                    wlast.rm = clast.rm;
                end
                ax = subplot(r, c, c*(i-1) + j-1,'Parent', m.parent);
                cc = m.crossed{m.wini,i,j};%xcorr2(w1.rm,w2.rm);
                imagesc(ax, cc); title(ax, sprintf('c%d x c%d',w1.ind,w2.ind)); %upside down?
                axis(ax,'equal'); axis(ax,'tight'); set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
                xlabel(ax, 'xcorr');
                colormap(ax,'jet');
                hold(ax,'on'); plot(ax, size(cc,2)/2,size(cc,1)/2,'md','MarkerFaceColor','w','MarkerSize',7)
            end
            %plot rate mat
            ax = subplot(r,c,i*r,'Parent', m.parent); imagesc(ax, w1.rm);
            xlabel(ax,'RM');title(ax,sprintf('c%d',w1.ind));
            axis(ax,'equal'); axis(ax,'tight'); set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
            myColorMap = jet(256); myColorMap(1,:) = 1; colormap(ax, myColorMap);
            %colormap(ax,'jet');
        end
        %plot last rate mat ?
        ax = subplot(r,c,r*r,'Parent', m.parent); imagesc(ax, wlast.rm);
        xlabel(ax,'RM');title(ax, sprintf('c%d',wlast.ind));
        axis(ax,'equal'); axis(ax,'tight'); set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
        myColorMap = jet(256); myColorMap(1,:) = 1; colormap(ax, myColorMap);
        % SCATTER of gridscore
        ax = subplot(r,c,r*(r-1)+1,'Parent', m.parent,'Position',[0.05 0.09 0.35 0.45]); 
        plot(ax, 1:length(m.g),[a; b],'o'); legend(ax,{'before','muscimol'});
        xlabel(ax,' cell');ylabel(ax,'gridscore');title(ax, sprintf('%s','group gridscore after muscimol'));
        axis(ax,'equal'); ylim(ax,[-2 2]);  %set(ax,'ydir','normal','xticklabel',[],'yticklabel',[]);
    else
        %empty group
    end
    
end

