function plotGroup(args, ind, splitPages)
    %groups = findSimultaneouslyRecordedCells(args.cells);
    %k = fieldnames(groups);
    midmax = 0;
    for i = 1:length(args.groups)
        g = args.groups{i};
        g = sortByMiddle(g);
        midmax = max(midmax ,length(g(1).middle));
    end
    midmax = 6;
    if ~isfield(args,'parent')
        args.parent = figure('Position', [0, 0, 2000, 1000]);
        set(gca,'LooseInset', get(gca,'TightInset'));
    end
    if i == -1
        for i = 1:length(args.groups)
            fprintf('%s %d %d\n', 'i-', i, length(args.groups{i}));
            plotGroupPrivate(args.parent, args.groups{i}, args.grid_thresh, i, midmax, splitPages);
        %set(suptitle(ax, titl),'Interpreter', 'none'); %PUT SUP TITLE AFTER ALL SUBPLOT COMMANDS << test
        end
    else
        %test(args.parent);
        plotGroupPrivate(args.parent, args.groups{ind}, args.grid_thresh, ind, midmax, splitPages);        
        %set(suptitle(ax, titl),'Interpreter', 'none'); %PUT SUP TITLE AFTER ALL SUBPLOT COMMANDS << test
    end
end

function [ax ai] = axi(parent, ai, ncol, naxes)
    nRow = ceil( naxes / ncol ) ;
    rowH = 0.8 / nRow ;  colW = 0.8 / ncol ; %width of plot
           %offset
    colX = 0.05 + linspace( 0, 0.9, ncol+1 ) ; colX = colX(1:end-1) ;
    rowY = 0.05 + linspace( 0.9, 0, nRow+1 ) ; rowY = rowY(2:end);
    rowId = ceil( ai / ncol );
    colId = ai - (rowId - 1) * ncol;
    ax = axes('Parent',parent, 'Position', [colX(colId), rowY(rowId), colW, rowH] ) ;
    ai = ai + 1; 
end

function titl = plotGroupPrivate(parent, g, thresh, iter, midmax, splitPages)
    rmt = 1.5; % Rate Map Threshold Hz
    debug = '';
    good = [g(1)]; bad = [g(1)]; %find low gridscore cells
    for j = 1:length(g) %put low grid scores at end
        if g(j).before.gridscore > thresh
            good(end+1) = g(j);
        else
            bad(end+1) = g(j);
        end
    end
    g = good(2:end); %[good(2:end) bad(2:end)]; Only show good sells
    if ~isempty(g)
        %g = sortByMiddle(g); %align times of cells
        %g = sortByEllipse(g); %display by ellipse size;
    end
    ncol = 3 * (1+length(g(1).middle)+1) + 1;
    naxes = length(g) * ncol;

    %MAIN LOOP
    ai = 1;
    for i = 1:length(g)
        r1 = g(i);
        mcross = {};
        gc = '';
        if r1.before.gridscore < thresh
            gc = '*';
        end
        %Trajectory
        %ax = subplot(r, c, c*(z-1) + 1,'Parent', parent);
        [ax ai] = axi(parent, ai, ncol, naxes);
        plot(ax, r1.before.px, (r1.before.py));
        hold(ax, 'on');
        scatter(ax, r1.before.sx, (r1.before.sy), '.'),...
            xlim(ax, [0 100]), ylim(ax, [0 100]);
        title(ax, sprintf('%sc%d:tet%d: before',gc, r1.ind,r1.tet),'fontweight','bold');
        %plot middle trajectory
        for i = 1: min(length(r1.middle), midmax)
            if  goodCell(r1.middle{i}, rmt)
                %ax = subplot(r, c, c*(z-1) + 1 + ind, 'Parent', parent);
                [ax ai] = axi(parent, ai, ncol, naxes);
                plot(ax, r1.middle{i}.px, (r1.middle{i}.py));
                hold(ax, 'on');
                scatter(ax, r1.middle{i}.sx, (r1.middle{i}.sy), '.');
                xlim(ax, [0 100]); ylim(ax, [0 100]);
                title(ax, sprintf('%.0f-%.0f min', r1.middle{i}.pt(1)/60,...
                    r1.middle{i}.pt(length(r1.middle{i}.pt))/60));
            else
                ai = ai+1;
            end
            
        end                           %after
        if r1.after.exists == 1
            %ax = subplot('Parent', parent, r, c, c*(z-1) + 1 + m + 1);
            [ax ai] = axi(parent, ai, ncol, naxes);
            plot(ax, r1.after.px, (r1.after.py));
            hold(ax, 'on');
            scatter(ax, r1.after.sx, (r1.after.sy), '.'),...
                xlim(ax, [0 100]), ylim(ax, [0 100]);
            title(ax, sprintf('%s', 'after'));
        end
        
        %plot middle rms
        %ax = subplot(r, c, c*(z-1) + 1, 'Parent', parent);
        [ax ai] = axi(parent, ai, ncol, naxes);
        imagesc(ax, r1.before.rm); 
        title(ax, sprintf('gs%.1f mxHz%0.f',round(r1.before.gridscore,1), r1.before.max_r));
        mcross{1} = r1.before.rm; 
        for i = 1: min(length(r1.middle), midmax)
            if  goodCell(r1.middle{i}, rmt)
                %ax = subplot(r, c, c*(z-1) + 1 + ind, 'Parent', parent);
                [ax ai] = axi(parent, ai, ncol, naxes);
                imagesc(ax, r1.middle{i}.rm);
                title(ax, sprintf('gs%.1f mxHz%0.f',round(r1.middle{i}.gridscore,1), r1.middle{i}.max_r));
                mcross{end + 1} = r1.middle{i}.rm;
            else
                ai = ai+1;
            end
            
        end                    
        if r1.after.exists == 1
            [ax ai] = axi(parent, ai, ncol, naxes);
            imagesc(ax, r1.after.rm);
            title(ax, sprintf('gs%.1f mxHz%0.f',round(r1.after.gridscore,1), r1.after.max_r));
            mcross{end + 1} = r1.after.rm;
        end
        %%%plot ac
        [ax ai] = axi(parent, ai, ncol, naxes);
        plotAcorrModule(r1.before, ax, '');
        %plot middle ac
        for i = 1: min(length(r1.middle), midmax)
            if goodCell(r1.middle{i}, rmt)
                %ax = subplot('Parent', parent, r, c, c*(z-1) + 1 + ind);
                [ax ai] = axi(parent, ai, ncol, naxes);
                plotAcorrModule(r1.middle{i}, ax, '');
            else
                ai = ai+1;
            end
        end                           %after
        if r1.after.exists == 1
            [ax ai] = axi(parent, ai, ncol, naxes);
            plotAcorrModule(r1.after, ax, ''); %plot after mid
            
        end
        %mcross
        %ax = subplot(r, c, c*(z-1) + 2 + 2*midmax + 3,'Parent', parent);
        mc = zeros(length(mcross));
        mc = allCorr(mcross,mcross); %ADD BACK 
        mc =  padarray(mc,midmax+2-size(mc),-1,'post');
        mc = flip(mc); 
        [ax ai] = axi(parent, ai, ncol, naxes);
        imagesc(ax, mc);caxis(ax, [-1 1]);
        title(ax, sprintf('[%.1f %.1f]', min(mc(:)), max(mc(:))));
        %time correlation
        %           ax = subplot('Parent', parent, r, c, c*(z-1) + 3 + m);
        %           t2 = g(j).before.st;
        %           [p, c] = time_correlation(t1,t2);
        %MAKE AXIS NICE
        ax = findobj(parent,'Type','Axes');
        for i=1:length(ax)
            set(ax(i),'LooseInset', get(ax(i),'TightInset'));
            set(ax(i),'FontSize',8);
            axis(ax(i),'equal')
            axis(ax(i),'tight');
            axis(ax(i),'off')
            colormap(ax(i),'jet')
            set(ax(i),'ydir','normal');
            %axis off; axis equal;
        end
        
    end
    %ADD TITLE
    ax = axes('Parent',parent, 'Position', [.1, .95, 1, .02],'FontSize',9, 'fontweight','bold'); 
    axis(ax,'off');
    txt = 'Muscimol session [ Before(t-1)  Mids(t_0-t_n)  After(t+1) ]: Trajectory, Rate Maps,  Autocorrelations,  Cross correlation of all sessions';
    text(ax, 0,0, txt);
    
    
end

function g = sortByMiddle(g)
    fs = fieldnames(g);
    ie = find(strcmp(fs, 'middle'));
    a = struct2cell(g);
    s = size(a);
    a = reshape(a, s(1), []);
    a = a';
    for i = 1:length(a(:,1))
        a(i,length(fs)+1) = {length(a{i,ie})}; %middle row
    end
    a = sortrows(a, -length(fs)+1);
    a = a(:,1:length(fs));
    a = reshape(a', s);
    g = cell2struct(a, fs, 1);
end


function g = sortByEllipse(g)
    fs = fieldnames(g);
    ie = find(contains(fs, 'before'));
    a = struct2cell(g);
    s = size(a);
    a = reshape(a, s(1), []);
    a = a';
    for i = 1:length(a(:,1))
        %sort by ellipse size of *before*
        a(i,length(fs)+1) = {a{i,ie}.('module').('major_ax') * a{i,ie}.('module').('minor_ax')*pi};
        if isnan(a{i,ie}.('module').('major_ax'))  || isnan(a{i,ie}.('module').('minor_ax'))
            a(i,length(fs)+1) = {10000}; %make unknown max
        end
    end
    a = sortrows(a, length(fs)+1); %-11 for backwards
    a = a(:,1:length(fs));
    a = reshape(a', s);
    g = cell2struct(a, fs, 1);
end

function b=goodCell(c, rmt)
    b= c.exists && ...
        not(isnan(c.gridscore) || ...
        c.gridscore == -2 || ...
        c.max_r < rmt     || ...
        c.max_r == 50        ...
        );
end

function plotAcorrModule(c, ax, s)%m,n,l
    %ax = subplot(sub(1), sub(2), sub(3),'Parent', parent); %(m,n,l);
    imagesc(ax, c.ac');  hold(ax, 'on');%xlim(ax, xlim); ylim(ax, ylim);%axis('xy'); axis ij; axis equal; axis off;
    x=xlim(ax); y=ylim(ax);
    title(ax, sprintf('%sspk%d',s, length(c.sx)));
    lw = 0.5;
    major = c.module.major_ax; minor = c.module.minor_ax; hex_peaks = c.module.hex_peaks;
    phi = c.module.angle; %+ pi/2;
    co = 'k';
    %center point
    if hex_peaks ~= -1 %CHANGE TO EXISTS
        x0=hex_peaks(7,1); y0=hex_peaks(7,2);
        beta = phi; sinbeta = sin(beta); cosbeta = cos(beta);
        alpha =0: pi/100:2*pi; sinalpha = sin(alpha);cosalpha = cos(alpha);
        x1 = x0 + (major * cosalpha * cosbeta - minor * sinalpha * sinbeta);
        y1 = y0 + (major * cosalpha * sinbeta + minor * sinalpha * cosbeta);
        x1(x1<x(1)) = x(1); x1(x1>x(2)) = x(2);y1(y1<y(1)) = y(1); y1(y1>y(2)) = y(2);
        plot(ax, x1,y1,co,'LineWidth',lw); 
        plot(ax, hex_peaks(:,1),hex_peaks(:,2),'ok','LineWidth',lw); 
        xMajor1 = x0 + major * cos(phi); xMajor2 = x0 - major * cos(phi);
        yMajor1 = y0 + major * sin(phi); yMajor2 = y0 - major * sin(phi);
        p1=xMajor1:(xMajor2-xMajor1)/10:xMajor2; p2=yMajor1:(yMajor2-yMajor1)/10:yMajor2;
        if ~isempty(p1) && ~isempty(p2)
            p1(p1<x(1)) = x(1); p1(p1>x(2)) = x(2);p2(p2<y(1)) = y(1); p2(p2>y(2)) = y(2);
            plot(ax, p1,p2, co,'LineWidth',lw); 
        end
        xMinor1 = x0 + minor * cos(phi+pi/2); xMinor2 = x0 - minor * cos(phi+pi/2);
        yMinor1 = y0 + minor * sin(phi+pi/2); yMinor2 = y0 - minor * sin(phi+pi/2);
        p11=xMinor1:(xMinor2-xMinor1)/10:xMinor2; p21=yMinor1:(yMinor2-yMinor1)/10:yMinor2;
        if ~isempty(p11) && ~isempty(p21)
            p11(p11<x(1)) = x(1); p11(p11>x(2)) = x(2);p21(p21<y(1)) = y(1); p21(p21>y(2)) = y(2);
            plot(ax, p11,p21, co,'LineWidth',lw);
        end
    end
    xlim(ax, x); ylim(ax, y);
end