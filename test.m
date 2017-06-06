
function test(m, z)
    %set(gui.status, 'String', sprintf('group(%d) session(%s)',gui.m.gid, t));
    
    train = [];
    c = [];
    for i = 1:length(m.g)
        if strcmp(m.sesh,'before')
            c = m.g{i}.before;
            edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
        elseif strcmp(m.sesh,'after')
            c = m.g{i}.after;
            edges = [-Inf, mean([c.pt(2:end)'; c.pt(1:end-1)']), +Inf];
        else
            c = m.g{i}.middle{m.sesh};
            edges = [-Inf, mean([c.pt(2:end); c.pt(1:end-1)]), +Inf];
        end
        %if isnumeric(m.sesh)
        I = discretize(c.st, edges);
        train(i,:) = logical(zeros(1, length(c.pt)));
        train(i,I) = 1;
    end
    train = logical(train);
    ax = axes('Parent',m.parent);
    %figure; ax = gca; %DELETE
    xlim(ax,[0 120]); ylim(ax,[0 120])
    colors = jet(length(m.g));
    %DRAW LOOP
    hold(ax, 'on');
    axis(ax,'off');
    colormap(ax,colors);
    i = 1;
    tic
    while i <= z
        if m.stop
            m.stop = false;
            cla(ax);
            break
        end
        while(m.pause)
            pause(1/5);
        end
        %trajectory
        %plot(ax,c.px(i),c.py(i),'g.');
        %delete(ax);
        %ax = axes('Parent',m.parent);
        for j = 1:length(m.g)
        plot(ax,c.px(1:i),c.py(1:i),'g.');
        end
        xlim(ax, [0 120]);ylim(ax, [0 120]);
        %axis(ax,'off');
        %plot cells
        hold(ax, 'on');
        tx = c.px(1:i); ty = c.py(1:i);
        %tx = repmat(c.px(1:i),length(m.g),1);ty = repmat(c.py(1:i),length(m.g),1);
        %plot(ax, double(tx(train(:,1:i))), double(ty(train(:,1:i)))+0.3*j,'o');%'Color', colors);
        %
        for j = 1:length(m.g)
            %if train(j,i) && m.cells(j)
            if m.cells(j)
                %text(ax, double(c.px(i)), double(c.py(i)), num2str(m.g{j}.ind),...
                    %'fontsize',14,'fontweight','bold','Color', colors(j,:));
                %t = logical(zeros(1,length(c.pt)));t(1:i) = train(j,1:i);%t = logical(t);
                %text(ax, double(c.px(t)), double(c.py(t)), num2str(m.g{j}.ind),...
                    %'fontsize',14,'fontweight','bold','Color', colors(j,:));
            %text(ax, double(tx(train(j,1:i))), double(ty(train(j,1:i))),...
                %num2str(m.g{j}.ind),'fontsize',14,'fontweight','bold','Color', colors(j,:));
            plot(ax, double(tx(train(j,1:i))), double(ty(train(j,1:i)))+0.3*j,'o','Color', colors(j,:));
            end
        end
        %}
        text(ax, 110,110,num2str(i),'fontsize',14,'color','k');
        text(ax, 110,110,num2str(i),'fontsize',14,'color','w');
        hold(ax, 'off');
        drawnow();
        %delete(ax);
        pause(1/m.speed);
        i = i + 1;
    end
    toc
end%run()

function test1(c, z)
    i = 0;
    figure;
    hold off
    tic
    while i < z
        i = i+1;
        plot(c.px(1:i),c.py(1:i),'k.');
        %plot(c{1}.before.px(i),c{1}.before.py(i),'k.');
        xlim([0 120]);ylim([0 120]);
        text(110,110,num2str(i),'fontsize',14,'color','k');
        text(110,110,num2str(i),'fontsize',14,'color','w');
        drawnow();
        pause(1/20000);
    end; toc
end


function testGroup(parent)
    
    % - Define dummy data: 11 time series.
    t       = 0 : 0.1 : 10 ;
    data    = 2 * repmat( sin(t).', 1,11 ) + rand( length(t), 11 ) ;
    nSeries = size( data, 2 ) ;
    % - Build figure.
    %figure() ;  clf ;
    %set( gcf, 'Color', 'White', 'Unit', 'Normalized','Position', [0.1,0.1,0.6,0.6] ) ;
    % - Compute #rows/cols, dimensions, and positions of lower-left corners.
    %%
    nCol = 4 ;  nRow = ceil( nSeries / nCol ) ;
    rowH = 0.58 / nRow ;  colW = 0.7 / nCol ;
    colX = 0.06 + linspace( 0, 0.96, nCol+1 ) ;  colX = colX(1:end-1) ;
    rowY = 0.1 + linspace( 0.9, 0, nRow+1 ) ;  rowY = rowY(2:end) ;
    % - Build subplots axes and plot data.
    for dId = 1 : nSeries
        rowId = ceil( dId / nCol ) ;
        colId = dId - (rowId - 1) * nCol ;
        ax = axes('Parent',parent, 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
        plot(ax, t, data(:,dId), 'b' ) ;
        grid on ;
        xlabel(ax, '\theta(t) [rad]' ) ;  ylabel(ax, 'Anomaly [m]' ) ;
        title(ax, sprintf( 'Time series %d', dId )) ;
    end
    % - Build title axes and title.
    ax = axes('Parent',parent, 'Position', [0, 0.95, 1, 0.05] ) ;
    set(ax, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
    text(ax, 0.5, 0, 'My Nice Title', 'FontSize', 14', 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ) ;
end