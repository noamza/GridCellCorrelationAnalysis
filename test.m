
function test(m, z)
    fn = sprintf('C:\\Noam\\Data\\muscimol\\noam\\cells_%dmin_b_midscorrect.mat',10);    
    fprintf('loading %s ',fn); %ascii 48
    tic; cells = load(fn); cells = cells.cells; toc;
    
    c = cells{1};  c = c.before;
    
    figure; colormap jet; s  = 10
    ac = c.ac;
    subplot(2,1,1) 
    imagesc(ac);hold on;%%%%
    [zmax,imax,zmin,imin]= Extrema2(ac);
    [i,j]=ind2sub(size(ac),imax);
    plot(j,i,'wo','markersize',s);
    pks = FastPeakFind(ac);
    plot(pks(1:2:end),pks(2:2:end),'w+','markersize',s);
    ac(isnan(ac)) = 0;
    [k,l] = find(imregionalmax(ac));
    plot(l,k,'wx','markersize',s); %%%%%
    ac = xcorr2(a.rm,a.rm);
    subplot(2,1,2)  
    imagesc(ac);hold on;%%%%
    [zmax,imax,zmin,imin]= Extrema2(ac);
    [i,j]=ind2sub(size(ac),imax);
    plot(j,i,'wo','markersize',s);
    pks = FastPeakFind(ac);
    plot(pks(1:2:end),pks(2:2:end),'w+','markersize',s);
    ac(isnan(ac)) = 0;
    [k,l] = find(imregionalmax(ac));
    plot(l,k,'wx','markersize',s); %%%%%

        
        num_bins = 50;
        rate_map_time =   zeros(num_bins, num_bins);
        rate_map_spikes = zeros(num_bins, num_bins);
        for i = 1:length(c.pt)
            col =            ceil(c.px(i)/max(c.px) * num_bins);
            row =            ceil(c.py(i)/max(c.py) * num_bins); %flips y
            rate_map_time(row, col) = rate_map_time(row, col) + 1;%0.02;
        end
        a = []; hold on;
        for i = 1:length(c.st)
            col =            ceil(c.px(c.si(i))/max(c.px) * num_bins);
            row =            ceil(c.py(c.si(i))/max(c.py) * num_bins);
            rate_map_spikes(row, col) = rate_map_spikes(row, col) + 1;
            % because some spikes are in other positions by interpolation
            if(rate_map_time(row, col) < rate_map_spikes(row, col))
                disp('here');
                a(end+1,:) = [col row]; 
                 %rate_map_time(row, col) = rate_map_time(row, col) + 0.02;
            end
            %rate_map_time(row, col) = rate_map_time(row, col) + 1; %% IS THIS VALID? double counting???
        end
        rate_map = (rate_map_spikes ./ rate_map_time);
        rate_map(isnan(rate_map)) = 0;
        max(rate_map(:));
        rate_map = rate_map/max(rate_map(:)); %normalize
    
    

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