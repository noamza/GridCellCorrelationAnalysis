%%%%%%%%%
% N Z A % 
% grid cells: 228
% pairs: 749
%%% DATASTATS() !!
%%%%%%%%%

%{
    ****Questions***
- Nan errors in gridscore, etc
- Method for correlations that don't fit

% differences: mid length
               after existence
               after lengths
               area
               number of middle fields
               no overlapping tet/ch combo (up to 9?)
               g12048_030708 * 54956 31508 54956 *
               no change in before length
%}

function musciomol()
    
    params.dir_load =...
        'C:\Noam\Data\muscimol\DB_musc_MEC\';
    params.dir_save =....
        'C:\Noam\Data\muscimol\noam_output\';
    
%    tic %1800 sec   %LOAD FROM FILES
%     
     files = dir(strcat(params.dir_load,'DB*.mat'));
     for i=1:2%length(files) 
         data{i} = load(strcat(params.dir_load,files(i).name));
         data{i}.db.ind = i;
         cells{i} = process(data{i}.db);
         fprintf('%.1f \n',(i*100/length(files)));
     end
     %save(strcat('C:\Noam\Data\muscimol\noam\cells_10min.mat',''),'cells');
     %toc
     %}
    tic; cells = load('C:\Noam\Data\muscimol\noam\cells_10min.mat'); cells = cells.cells; toc;
    groups = find_simultaneously_recorded_cells(cells);
    k = fieldnames(groups);
    
    for i = 1:length(k);
       g = groups.(k{i});
       t = g(1); t= length(t.after.trajectory.px);
       %t.a; length(t.after.trajectory.px) .exists length(t.middle)
       fprintf('%s * %d ',k{i} ,t);
       %fprintf('%s %d', k{i}, length(g));
       for j = 1:length(g) %in group
            tt = g(j); tt = length(tt.after.trajectory.px);
            %.a length(tt.after.trajectory.px) .exists length(tt.middle)
            if t ~= tt               %not(strcmp(t,tt))
                fprintf('%d ', tt);
                t = tt;
            end; %fprintf('%d %d\n', g(j).tet, g(j).cel);
       end; disp('*'); end;
   
   
   
    %pairs
    tic
    pairs = find_pairs(cells);%find_pairs_of_cells_mosimol(data);
    d = 47; %65 23 86 235
    for i = 48:49 %47:length(pairs)
        r1i = pairs(i,1); r2i = pairs(i,2);
        fprintf('%.1f of %d ',i*100/length(pairs),length(pairs));
        r1 = cells{r1i};%process(data{r1i}.db);
        r2 = cells{r2i};%process(data{r2i}.db);
        plot_pair(r1, r2, i);
    end
    toc
    disp('muscimol done');
    %{
    - Gilad's paper
    - to what extent do cells disappear together
    - max simul recorded
    - module
    %}
    
    
    
    
end

function plot_group(groups)
figure('Position', [100, 100, 1200, 900]);
set(gca,'LooseInset',get(gca,'TightInset'));
colormap jet;
rmt = 0.1; % Rate Map Threshold
keys = fieldnames(groups);

for k = 1:length(keys)
    key = keys{k};
    g = groups.(key);
    set(suptitle(sprintf('%s',key)));    
    r = length(g); %Rows for figure (trajectory ratemap autoc xcross axcross)
    n = 99;
    
    %%%plot rm
    subplot(r, n, 1 + n*2);
    imagesc(r1.before.rm), title(sprintf('c1: %.0fHz', r1.before.max_r)); axis off; axis equal;
    subplot(r, n, 1 + n*3);
    imagesc(r2.before.rm), title(sprintf('c2: %.0fHz', r2.before.max_r)); axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*2);
        imagesc(r1.middle{i}.rm);title(sprintf('%.0fHz', r1.middle{i}.max_r));axis off; axis equal;
        subplot(r, n, i+1 + n*3);
        imagesc(r2.middle{i}.rm);title(sprintf('%.0fHz', r2.middle{i}.max_r));axis off; axis equal;
    end                           %after
    if after
        subplot(r, n, n*3-1);
        imagesc(r1.after.rm);title(sprintf('%.0fHz', r1.after.max_r));axis off; axis equal;
        subplot(r, n, n*4-1);
        imagesc(r2.after.rm);title(sprintf('%.0fHz', r2.after.max_r));axis off; axis equal;
    end
    % mcross
    t1{1} = r1.before.rm; t2{1} = r2.before.rm;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.rm; t2{i+1} = mid2{i}.rm;
    end
    if after
        t1{m+2} = r1.after.rm; t2{m+2} = r2.after.rm;
    end
    subplot(r, n, n*3);
    imagesc(allCorr(t1,t1));title('cross per');axis off; axis equal;
    subplot(r, n, n*4);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;
    
    %%%plot ac
    subplot(r, n, 1 + n*4);
    imagesc(r1.before.ac), title(sprintf('c1: %.2f', r1.before.gridscore));
    axis off; axis equal;
    subplot(r, n, 1 + n*5);
    imagesc(r2.before.ac), title(sprintf('c2: %.2f', r2.before.gridscore));
    axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*4);
        imagesc(r1.middle{i}.ac); title(sprintf('%.2f', r1.middle{i}.gridscore)); axis off; axis equal;
        subplot(r, n, i+1 + n*5);
        imagesc(r2.middle{i}.ac); title(sprintf('%.2f', r2.middle{i}.gridscore)); axis off; axis equal;
    end
    if after
        subplot(r, n, n*5-1);
        imagesc(r1.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r1.after.gridscore));
        subplot(r, n, n*6-1);
        imagesc(r2.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r2.after.gridscore));
    end
    % mcross
    t1 = {}; t2 = {}; t1{1} = r1.before.ac; t2{1} = r2.before.ac;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.ac; t2{i+1} = mid2{i}.ac;
    end
    if after
        t1{m+2} = r1.after.ac; t2{m+2} = r2.after.ac;
    end
    subplot(r, n, n*5);
    imagesc(allCorr(t1,t1)); title('cross per');axis off; axis equal;
    subplot(r, n, n*6);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;
    
    %%% cross corr
    subplot(r, n, 1 + n*6);
    cc = Cross_Correlation(r1.before.rm, r2.before.rm);
    imagesc(cc); title(sprintf('Xcorr: %s',''));axis off; axis equal;
    %acc
    subplot(r, n, 1 + n*7);
    acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
    imagesc(acc);
    title(sprintf('AcXcorr: %.2f', gridscore(acc, -1))); axis off; axis equal;
    for i = 1: m;
        subplot(r, n, i+1 + n*6);
        cc = Cross_Correlation(r1.middle{i}.rm, r2.middle{i}.rm);
        imagesc(cc);axis off; axis equal;
        %acc
        subplot(r, n, i+1 + n*7);
        acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
        imagesc(acc);axis off; axis equal;title(sprintf('%.2f', gridscore(acc, -1)));
    end
    
    if after
        subplot(r, n, n*7-1);
        cc = Cross_Correlation(r1.after.rm, r2.after.rm); imagesc(cc);axis off; axis equal;
        %acc
        subplot(r, n, n*8-1); acc = Cross_Correlation(cc, cc);
        getSubmatrixFromCenter(acc, cc); imagesc(acc);
        title(sprintf('%.2f', gridscore(acc, -1)));axis off; axis equal;
    end
    
    
    
    
    sdir = 'C:\Noam\Output\muscimol\';
    filename = sprintf('%s%d_Rat_%s_date_%s_C%d_t%d_c%d_C%d_t%d_c%d.png',sdir, id, r1.id, r1.date, r1.ind, r1.tet, r1.cel, r2.ind, r2.tet, r2.cel);
    disp(filename);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(filename, '-dpng','-r0');
    close;
end
end
    
    
function groups = find_simultaneously_recorded_cells(cells)
    groups = [];
    gi = 1;
    for i=1:length(cells) 
        a = cells{i};
        %if not(a.max_r < rate_map_thresh || a.max_r == 50)
        key = sprintf('g_%s_%s', a.id, a.date); 
        if isfield(groups, key)
            t = groups.(key);
            t = [t a];
            groups.(key) = t;
        else
            groups = setfield(groups, key, a);
        end
        %    fprintf('%d: bad cell %f rmthr %f', i, a.max_r, rate_map_thresh);
    end
    
    %validate
    keys = fieldnames(groups);
    for i = 1:length(keys)
        k = keys{i}; 
        group = groups.(k); 
        t = group(1); 
        t = length(t.before.trajectory.px);
        for j = 1:length(group) %CHECK TIMES !!
            tt = group(j); %cells{group(j)};
            tt = length(tt.before.trajectory.px);
            if t ~= tt
                fprintf('%f %f whats up with this group b.length %s?\n',t, tt, k);
            end
        end
    end
    
    fprintf('grouped\n');
end

function plot_pair(r1, r2, id)
    figure('Position', [100, 100, 1200, 900]);
    set(gca,'LooseInset',get(gca,'TightInset'));
    colormap jet; bad = '';
    rmt = 0.1; % Rate Map Threshold
    
    %n = min(length(r2.middle), length(r1.middle));
    if length(r2.middle) ~= length(r1.middle)
        %fprintf('%d: error plotting, middle sessions different lengths', id);%bad = '*';
    end
    
    set(suptitle(sprintf('%sP%d(i%d,i%d): rat %s %s C1(t%dc%d %s %s) C2(t%dc%d %s %s)',...
    bad,id,r1.ind,r2.ind,r1.id,r1.date,r1.tet,r1.cel,r1.a,r1.type,r2.tet,r2.cel,r2.a,r2.type)),...
    'Interpreter', 'none');
    
    %line up bin times
    i = 1; i1 = 1; i2 = 1; mid1 = {}; mid2 = {};
    while i1 <= length(r1.middle) && i2 <= length(r2.middle)
        if r1.middle{i1}.pt(1) == r2.middle{i2}.pt(1) %Will only show traj,ratemap/ac for good ratemaps
            if not(r1.middle{i1}.max_r < rmt || r1.middle{i1}.max_r == 50 ...
            || r2.middle{i2}.max_r < rmt || r2.middle{i2}.max_r == 50)
                mid1{i} = r1.middle{i1};
                mid2{i} = r2.middle{i2};
                i = i + 1;
            end
            i1 = i1 + 1;
            i2 = i2 + 1;
        elseif r1.middle{i1}.pt(1) < r2.middle{i2}.pt(1)
            i1 = i1 + 1;
        else
            i2 = i2 + 1;
        end
    end
    r1.middle = mid1;
    r2.middle = mid2;
    m = length(mid1);
    n = m +1+1; %how many cols in figure, +1 before +1 mcross
    after = 0; %display after;
    if not(r1.after.max_r < rmt || r1.after.max_r == 50 ...
         ||r2.after.max_r < rmt || r2.after.max_r == 50)
        n = n + 1;
        after = 1;
    end
    r = 8; %Rows for figure (trajectory ratemap autoc xcross axcross)
    
    
    %%%plot trajectory
    subplot(r, n, 1 + n*0);
    plot(r1.before.trajectory.px, flip(r1.before.trajectory.py));
    hold on
    scatter(r1.before.trajectory.sx, flip(r1.before.trajectory.sy), '.'),...
        xlim([0 100]), ylim([0 100]), ...
        title(sprintf('c1: before %d',length(r1.before.trajectory.st))); 
    axis off; axis equal;
    subplot(r, n, 1 + n*1);
    plot(r2.before.trajectory.px, flip(r2.before.trajectory.py));
    hold on
    scatter(r2.before.trajectory.sx, flip(r2.before.trajectory.sy), '.'),...
        xlim([0 100]), ylim([0 100]), axis off; axis equal; 
    title(sprintf('c2: %d',length(r2.before.trajectory.st))); %TITLE
    for i = 1:m;
        %if length(r1.middle{i}.px) > 1 && length(r2.middle{i}.px) > 1 %Display if more than one spike
            subplot(r, n, i+1 + n*0);
            plot(r1.middle{i}.px, flip(r1.middle{i}.py));
            hold on
            scatter(r1.middle{i}.sx, flip(r1.middle{i}.sy), '.'); xlim([0 100]), ylim([0 100]);
                %title(sprintf('%.0f-%.0f(%d)', min(r1.middle{i}.pt)/60,... 
                        %max(r1.middle{i}.pt)/60, length(r1.middle{i}.st)));
            title(sprintf('%.0f-%.0f(%d)', r1.middle{i}.pt(1)/60,... 
                        r1.middle{i}.pt(length(r1.middle{i}.pt))/60, length(r1.middle{i}.st)));
            axis off; axis equal;
            subplot(r, n, (i+1) + n*1);
            plot(r2.middle{i}.px, flip(r2.middle{i}.py));
            hold on
            scatter(r2.middle{i}.sx, flip(r2.middle{i}.sy), '.'), xlim([0 100]), ylim([0 100]);
            %title(sprintf('%d',length(r2.middle{i}.st)));
            title(sprintf('%.0f-%.0f(%d)', r2.middle{i}.pt(1)/60,... 
                        r2.middle{i}.pt(length(r2.middle{i}.pt))/60, length(r2.middle{i}.st)));
            axis off; axis equal;
    end
    %if length(r1.after.trajectory.px) > 1 && length(r2.after.trajectory.px) > 1
    if after
        subplot(r, n, n*1-1);
        plot(r1.after.trajectory.px, flip(r1.after.trajectory.py));
        hold on
        scatter(r1.after.trajectory.sx, flip(r1.after.trajectory.sy), '.'),...
            xlim([0 100]), ylim([0 100]), axis off; axis equal;
        title(sprintf('after %d',length(r1.after.trajectory.st))); 
        subplot(r, n, n*2-1);
        plot(r2.after.trajectory.px, flip(r2.after.trajectory.py));
        hold on
        scatter(r2.after.trajectory.sx, flip(r2.after.trajectory.sy), '.'),...
            xlim([0 100]), ylim([0 100]); axis off; axis equal;
        title(sprintf('%d',length(r2.after.trajectory.st)));
    end
    
    %%%plot rm 
    subplot(r, n, 1 + n*2);
    imagesc(r1.before.rm), title(sprintf('c1: %.0fHz', r1.before.max_r)); axis off; axis equal;
    subplot(r, n, 1 + n*3);
    imagesc(r2.before.rm), title(sprintf('c2: %.0fHz', r2.before.max_r)); axis off; axis equal;
    for i = 1: m;
            subplot(r, n, i+1 + n*2);
            imagesc(r1.middle{i}.rm);title(sprintf('%.0fHz', r1.middle{i}.max_r));axis off; axis equal;
            subplot(r, n, i+1 + n*3);
            imagesc(r2.middle{i}.rm);title(sprintf('%.0fHz', r2.middle{i}.max_r));axis off; axis equal;
    end                           %after
    if after
        subplot(r, n, n*3-1);
        imagesc(r1.after.rm);title(sprintf('%.0fHz', r1.after.max_r));axis off; axis equal;
        subplot(r, n, n*4-1);
        imagesc(r2.after.rm);title(sprintf('%.0fHz', r2.after.max_r));axis off; axis equal;
    end
    % mcross
    t1{1} = r1.before.rm; t2{1} = r2.before.rm;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.rm; t2{i+1} = mid2{i}.rm;
    end
    if after
        t1{m+2} = r1.after.rm; t2{m+2} = r2.after.rm; 
    end
    subplot(r, n, n*3);
    imagesc(allCorr(t1,t1));title('cross per');axis off; axis equal;
    subplot(r, n, n*4);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;
    
    %%%plot ac 
    subplot(r, n, 1 + n*4);
    imagesc(r1.before.ac), title(sprintf('c1: %.2f', r1.before.gridscore));
    axis off; axis equal;
    subplot(r, n, 1 + n*5);
    imagesc(r2.before.ac), title(sprintf('c2: %.2f', r2.before.gridscore));
    axis off; axis equal;
    for i = 1: m;
            subplot(r, n, i+1 + n*4);
            imagesc(r1.middle{i}.ac); title(sprintf('%.2f', r1.middle{i}.gridscore)); axis off; axis equal;
            subplot(r, n, i+1 + n*5);
            imagesc(r2.middle{i}.ac); title(sprintf('%.2f', r2.middle{i}.gridscore)); axis off; axis equal;
    end
    if after
        subplot(r, n, n*5-1);
        imagesc(r1.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r1.after.gridscore));
        subplot(r, n, n*6-1);
        imagesc(r2.after.ac); axis off; axis equal;
        title(sprintf('%.2f', r2.after.gridscore));
    end
    % mcross
    t1 = {}; t2 = {}; t1{1} = r1.before.ac; t2{1} = r2.before.ac;
    for i = 1:length(mid1)
        t1{i+1} = mid1{i}.ac; t2{i+1} = mid2{i}.ac;
    end
    if after
        t1{m+2} = r1.after.ac; t2{m+2} = r2.after.ac; 
    end
    subplot(r, n, n*5);
    imagesc(allCorr(t1,t1)); title('cross per');axis off; axis equal;
    subplot(r, n, n*6);
    imagesc(allCorr(t2,t2));title('cross per');axis off; axis equal;
    
    %%% cross corr
    subplot(r, n, 1 + n*6);
    cc = Cross_Correlation(r1.before.rm, r2.before.rm);
    imagesc(cc); title(sprintf('Xcorr: %s',''));axis off; axis equal;
    %acc
    subplot(r, n, 1 + n*7);
    acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
    imagesc(acc);
    title(sprintf('AcXcorr: %.2f', gridscore(acc, -1))); axis off; axis equal; 
    for i = 1: m;
            subplot(r, n, i+1 + n*6);
            cc = Cross_Correlation(r1.middle{i}.rm, r2.middle{i}.rm);
            imagesc(cc);axis off; axis equal;
            %acc
            subplot(r, n, i+1 + n*7);
            acc = Cross_Correlation(cc, cc); acc = getSubmatrixFromCenter(acc, cc);
            imagesc(acc);axis off; axis equal;title(sprintf('%.2f', gridscore(acc, -1)));
    end
    
    if after
        subplot(r, n, n*7-1);
        cc = Cross_Correlation(r1.after.rm, r2.after.rm); imagesc(cc);axis off; axis equal;
        %acc
        subplot(r, n, n*8-1); acc = Cross_Correlation(cc, cc);
        getSubmatrixFromCenter(acc, cc); imagesc(acc); 
        title(sprintf('%.2f', gridscore(acc, -1)));axis off; axis equal;
    end
    
    
    
    
    sdir = 'C:\Noam\Output\muscimol\';
    filename = sprintf('%s%d_Rat_%s_date_%s_C%d_t%d_c%d_C%d_t%d_c%d.png',sdir, id, r1.id, r1.date, r1.ind, r1.tet, r1.cel, r2.ind, r2.tet, r2.cel);
    disp(filename);
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    print(filename, '-dpng','-r0');
    close;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end plot
function f = flip(x)
    m = max(x);
    f = m-x;
end

%not perfect, submatrix will be a dimension higher than B, if odd
%(could fix using length instead of round(B)
function S = getSubmatrixFromCenter(A, B) %A big matrix, B little matrix
    raca = round(size(A)/2);
    rbcb = round(size(B)/2);
    S = A((raca(1)-rbcb(1)+1:raca(1)+rbcb(1)),((raca(2)-rbcb(2)+1:raca(2)+rbcb(2))));  
end

function A = allCorr(a, b)
   A = zeros(length(a));
   for i = 1:length(a);
       for j = 1:length(b)
           t1 = a{i}(:); t2 = b{j}(:);
           t1(isnan(t1))=0; t2(isnan(t2))=0;
           if length(t1) > length(t2)
               %could use image stretch <<
               t2 = [t2; zeros(length(t1)-length(t2),1)];
           else
               t1 = [t1; zeros(length(t2)-length(t1),1)];
           end
               A(i,j) = corr(t1, t2);
           if length(a{i}(:)) == length(b{j}(:))
               %A(i,j) = corr(a{i}(:), b{j}(:));
           else
               %A(i,j) = 1;
               %fprintf('wrong lengths %d %d\n',length(a{i}(:)),length(b{j}(:)));
           end
       end
   end
end



function r = process(data)
    r.id   = data.rat;
    r.date = data.date;
    r.tet  = data.tetrode;
    r.cel  = data.cell;
    r.type = data.cell_type_obj_new;
    r.a = data.area;
    r.ind = data.ind;
    p = data.B(1).pos_data;
    s = data.B(1).spike_data;
    %before
    r.before.trajectory = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
    tr = r.before.trajectory;
    r.before.rm = Create_Rate_Map(tr.px, tr.py, tr.pt, tr.sx, tr.sy, tr.st);
    r.before.ac = Cross_Correlation(r.before.rm, r.before.rm);
    r.before.max_r = max(r.before.rm(:));
    r.before.gridscore = gridscore(r.before.ac, r.ind);
    
    %middle
    p = data.B(2).pos_data;
    s = data.B(2).spike_data;
    tr = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
    binLength = 10; %minutes,
    maxTime =  45;% minutes
    bins = bin_trajactory(tr, binLength*60, maxTime*60); %BIN LENGTH
    for i = 1:length(bins)
        t = bins{i};
        bins{i}.rm = Create_Rate_Map(t.px, t.py, t.pt, t.sx, t.sy, t.st);
        bins{i}.ac = Cross_Correlation(bins{i}.rm, bins{i}.rm);
        bins{i}.max_r = max(bins{i}.rm(:));
        %fprintf('%d %d\n', r.ind, i);
        bins{i}.gridscore = gridscore(bins{i}.ac, r.ind);
    end
    r.middle = bins;
    
    %after
    if length(data.B) == 3               %after exists
        p = data.B(3).pos_data;
        s = data.B(3).spike_data;
        r.after.trajectory = rat_trajectory(p.x1, p.x2, p.y1, p.y2, p.t, s.x, s.x2, s.y, s.y2, s.ts);
        tr = r.after.trajectory;
        r.after.rm = Create_Rate_Map(tr.px, tr.py, tr.pt, tr.sx, tr.sy, tr.st);
        r.after.ac = Cross_Correlation(r.after.rm, r.after.rm);
        r.after.max_r = max(r.after.rm(:));
        r.after.gridscore = gridscore(r.after.ac, r.ind);
        r.after.exists = 1; 
    else                                %after doesn't exist
        r.after.trajectory = rat_trajectory([0],[0],[0],[0],[0],[0],[0],[0],[0],[0]);
        tr = r.after.trajectory;
        r.after.rm = Create_Rate_Map(tr.px, tr.py, tr.pt, tr.sx, tr.sy, tr.st);
        r.after.ac = Cross_Correlation(r.after.rm, r.after.rm);
        r.after.max_r = max(r.after.rm(:));
        r.after.gridscore = gridscore(r.after.ac, r.ind);
        r.after.exists = 1; 
        r.after.exists = false;
    end
    
end

function score = gridscore(ac, ind)
    try
        score = GridnessRadius(ac, FindROuter(ac)); 
    catch ME
        fprintf('%d: error calculting grid score %s\n', ind, ME.message);
        score = -2;
    end
end

function test()
    sa=0, sb= 0;
    for i = 1:length(bins)
        fprintf('%d: min %.3f max %.3f \n',i, min(bins{i}.pt)/60, max(bins{i}.pt)/60);
        sa = sa + length(bins{i}.pt);
        sb = sb + length(bins{i}.st);
    end
    fprintf('%d  %d \n',sa, sb);
end

%bin size in seconds, max time
%tr = trajectory; bs = bin size in seconds; mt = max size
function bins = bin_trajactory(tr, bs, mt)%
     bins{1,1} = []; b = 1;
     ind = 0;
     t = tr.pt; st = tr.st;
     t0 = ceil(tr.pt(1)/(5*60))*5*60 - bs; %stars at 5min intervals %tr.pt(1);
     for i = 1: length(t)
        if t(i) <= mt
            %breaks up bins based on bin length or 15min interval cutoffs
            if t(i) < t0 + bs %&& mod(t(i), bs) %mod breaks up bins on even intervals, non 0 means continue
                ind = ind +1;
                bins{b}.px(ind) = tr.px(i);
                bins{b}.py(ind) = tr.py(i);
                bins{b}.pt(ind) = tr.pt(i);
            else
                ind = 1;
                t0 = t(i);
                b = b + 1;
                bins{b}.px(ind) = tr.px(i);
                bins{b}.py(ind) = tr.py(i);
                bins{b}.pt(ind) = tr.pt(i);
            end
        end
     end
     %spikes
     b = 1; ind = 0;
     for i = 1:length(st)
        if st(i) <= mt && st(i) <= max(bins{length(bins)}.pt) %cuts off spikes past max trajectory time
            while st(i) > max(bins{b}.pt) %skips past bins with no spikes
                b = b + 1;
                ind = 0;
            end
            ind = ind +1;
            bins{b}.sx(ind) = tr.sx(i);
            bins{b}.sy(ind) = tr.sy(i);
            bins{b}.st(ind) = tr.st(i); 
        end
     end
     %empty bins
     minSpikes = 10;
     b = 1;
     t = {};
     for i = 1:length(bins)
        if isfield(bins{i},'sx') == 1 && length(bins{i}.st) > minSpikes; %if isfield(bins{i},'sx') == 0
            %bins{i}.sx(1) = 0;
            %bins{i}.sy(1) = 0;
            %bins{i}.st(1) = min(bins{b}.pt); 
            t{b} = bins{i};
            b = b + 1;
        end
     end
     bins = t;
end

function c = rat_trajectory(px1, px2, py1, py2, pt, sx1, sx2, sy1, sy2, st)
    c.px =  mean([px1, px2], 2)';
    c.px = c.px - min(c.px) + 0.00001; %no zeros
    c.py =  mean([py1, py2], 2)';
    c.py = c.py - min(c.py) + 0.00001;
    c.pt = pt;
    c.sx =  mean([sx1, sx2], 2)';
    c.sx = c.sx - min(c.sx) + 0.00001;
    c.sy =  mean([sy2, sy2], 2)';
    c.sy =c.sy - min(c.sy) + 0.00001;
    c.st = st;  
        %{figure();plot(c.px, c.py);hold on; plot(c.sx, c.sy,'.');%}
end

function rate_mat = Create_Rate_Map(px, py, pt, sx, sy, st)
    posx = px; posy = py; post = pt;
    spkx = sx; spky = sy; spkt = st;
    parms.sigma = 3; % gaussiam smoothing factor
    parms.time_per_bin=0.02; %defult size of bin time sampeling (1/sampeling rate)(updated in Read_Examples_2 function)
    % Minimum radius used in the auto-correlogram when finding the best
    parms.bin_size = 3; % size of spacial bin (for create the rate map)
    parms.num_of_direction_bins=60; % for head-direction calculatino
    parms.max_lag=500; % max lag (in msec) for temporal autocorrelation
    parms.bin_num=50;

    max_x = ceil(max(posx));
    max_y = ceil(max(posy));
    min_x = min(floor(posx));
    min_y = min(floor(posy));

    % divid the environment into spatial bins
    axis_x = min_x:parms.bin_size:max_x;
    axis_y = min_y:parms.bin_size:max_y;

    time_mat = zeros(length(axis_y),length(axis_x));
    spike_mat = zeros(length(axis_y),length(axis_x));
    rate_mat = zeros(length(axis_y),length(axis_x));

    %create time mat (compute how much time the rat spends in each bin)
    % find in each moment(time_per_bin) what spatial bin the rat is at and add the time_per_bin to
    %
    for i = 1:length(post)
        if ~isnan(posx(i)) && ~isnan(posy(i))
            [min_val,x_ind] =  min(abs(posx(i) - axis_x));
            [min_val,y_ind] =  min(abs(posy(i) - axis_y));
            time_mat(y_ind,x_ind) = time_mat(y_ind,x_ind) + parms.time_per_bin;

        end
    end
    %create conut mat( count the num of spikes in each bin)
    for i = 1:length(spkt)
        [min_val,x_ind] =  min(abs(spkx(i)-axis_x));
        [min_val,y_ind] =  min(abs(spky(i)-axis_y));
        spike_mat(y_ind,x_ind)= spike_mat(y_ind,x_ind)+1;
    end
    % create rate mat
    rate_mat=spike_mat./time_mat;
    rate_mat(rate_mat==inf)=NaN;
    rate_mat_unsmoothed = rate_mat;
    %create window
    h=fspecial('gaussian',3*[parms.sigma,parms.sigma],parms.sigma);
    rate_mat = nanconv2(rate_mat,h);
    rate_mat(isnan(rate_mat))=0; %bad??????nza
    disp('');
end

function pairs_by_file_index = find_pairs(cells)
    pairs_by_file_index = [];
    grid_cells = 0;
    for i=1:length(cells) %a = load(strcat(params.dir_load,files(i).name));
        a = cells{i};
        %if (strcmp(a.cell_type_obj_new, 'GRID') || strcmp(a.cell_type_obj_new, 'CONJ')) % '_GRID' ?????
        if a.before.gridscore > 0.3
            for j=i+1:length(cells)
                b = cells{j};
                %adjust for missing third <<<<<<
                if (       b.before.gridscore > 0.3 ...
                        && strcmp(a.id, b.id) && strcmp(a.date, b.date)...
                    ... && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                        && strcmp(a.a, b.a)  ...  %area????
                        && ( a.tet ~= b.tet || a.cel ~= b.cel) ...
                    ... && length(a.B) == 3 &&  length(b.B) == 3 ... %ignore cells without before/after
                    )
                    %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                        %a.date, a.tetrode, a.cell, b.tetrode, b.cell);

                    grid_cells = grid_cells + 1;
                    pairs_by_file_index(grid_cells,1) = i; 
                    pairs_by_file_index(grid_cells,2) = j; 
                end
            end
        else
            %disp (a.cell_type_subj);
        end
    end
    fprintf('%d pairs found\n',grid_cells);
end

function pairs_by_file_index = find_pairs_of_cells_mosimol(data)
    pairs_by_file_index = [];
    grid_cells = 0;
    for i=1:length(data) %a = load(strcat(params.dir_load,files(i).name));
        a = data(i); a = a{1}.db;
        if (strcmp(a.cell_type_obj_new, 'GRID') || strcmp(a.cell_type_obj_new, 'CONJ')) % '_GRID' ?????      
            for j=i+1:length(data)
                b = data(j); b = b{1}.db;
                if length(a.B) ~= 3
                    disp('missing 3rd');
                end
                if ( strcmp(a.rat, b.rat) && strcmp(a.date, b.date)...
                        && (strcmp(a.cell_type_subj, 'GRID') || strcmp(a.cell_type_subj, 'CONJ'))  ...
                        && strcmp(a.area, b.area)  ...  %area????
                        && ( a.tetrode ~= b.tetrode || a.cell ~= b.cell) ...
                        && length(a.B) == 3 &&  length(b.B) == 3 ... %ignore cells without before/after
                    )
                    %fprintf('same rat:%s date:%s a:%d,%d b:%d,%d \n', a.rat, ...
                        %a.date, a.tetrode, a.cell, b.tetrode, b.cell);

                    grid_cells = grid_cells + 1;
                    pairs_by_file_index(grid_cells,1) = i; 
                    pairs_by_file_index(grid_cells,2) = j; 
                end
            end
        else
            %disp (a.cell_type_subj);
        end
    end
    fprintf('%d paris found\n',grid_cells);
end

function R_outer = FindROuter(acorr)
    if length(acorr) == 1 || isnan(max(acorr(:)));
        R_outer = -1;
        %return;   %EXIT
    end
    % calculate all the extrema points in the spatial autocorrelation 
    [zmax,imax,zmin,imin]= Extrema2(acorr);
    [i,j]=ind2sub(size(acorr),imax);
    %put all extrema points in dist
    dist(:,1)=j;
    dist(:,2)=i;
    n=length(i);
    %calculate the distance of all extrema to the central peak and put them in
    %column 3
    dist(1:n,3)=sqrt(  (i(1:n)-i(1)).^2 + (j(1:n)-j(1)).^2);
     % sort the hexonal peaks by distance to the centeral peak
     [score,ind]=sort(dist(:,3));
     dist=dist(ind,:);
     %zmax=zmax(ind);
     R=dist(2,3);
      count=1;
     i=2;
     hex_peaks(1,:,:)=dist(1,:,:);
    % finds the first 6 closest peaks to the central peak
    while count<7 && i<=length(dist)
        % calculate the min distance of point i from all other already chosen
        % points
        min_dist_peaks=min(sqrt(  (hex_peaks(1:size(hex_peaks,1),1)-dist(i,1)).^2 +...
         (hex_peaks(1:size(hex_peaks,1),2)-dist(i,2)) .^2));
        % point i needs to be on the right side (we choos only half the point cause its semetrical)
        % and the distance of point i from all other already chosen points
        % needs to be higher than R/2
        if dist(i,1)>=dist(1,1) && min_dist_peaks>(R/1.5)
            hex_peaks(count,:,:)=dist(i,:,:); 
            count=count+1; 
            hex_peaks(count,1)=dist(1,1)-dist(i,1)+dist(1,1);
            hex_peaks(count,2)=dist(1,2)-dist(i,2)+dist(1,2);

            count=count+1;
        end    
        i=i+1;
    end
    R_outer=max(hex_peaks(:,3))*1.15;     
end

function gridness2 = GridnessRadius(org_Acor ,R_outer)
    if length(org_Acor) == 1 || isnan(max(org_Acor(:)))
        gridness2 = -1;
        %return; EXIT
    end
    % compute distance from center
    [val,center_y]=max(max(org_Acor));
    [val,center_x]=max(max(org_Acor'));
    [Y,X] = ndgrid(1:1:size(org_Acor,1), 1:1:size(org_Acor,2));
    dist_from_center=sqrt((Y-center_y).^2+(X-center_x).^2);
    % making sure that outer radius of the ring (R_outer) is not bigger than the distance matrix
    R_outer=min([min(dist_from_center(1,:)),min(dist_from_center(:,1)),...
        min(dist_from_center(size(dist_from_center,1),:))...
    min(dist_from_center(:,size(dist_from_center,2))),R_outer]);
    % compute inner radius of the anulus (ring)
    R_inner=ceil(min(dist_from_center(org_Acor<0.1))); 
    %extract the original anulus (ring) from Acor
    org_Ring=org_Acor(dist_from_center<=R_outer & dist_from_center>=R_inner);
    % make sure that after rotation and interpulation the center will remain the maximum point.
        org_Acor(center_x,center_y)=10;
        for jj = 2:6         
           % rotate the auto-correlation 
           rot_Acor=imrotate(org_Acor,(jj-1)*30,'bicubic');                 
             % compute distance from new center
            [val,tmp_center_x]=max(max(rot_Acor));
            [val,tmp_center_y]=max(max(rot_Acor'));
            [Y,X] = ndgrid(1:1:size(rot_Acor,1), 1:1:size(rot_Acor,2));
            tmp_dist_from_center=sqrt((Y-tmp_center_y).^2+(X-tmp_center_x).^2);
            % extract the anulus(ring)
            rot_Ring=rot_Acor(tmp_dist_from_center<=R_outer & tmp_dist_from_center>=R_inner);
            if length(rot_Ring)~=length(org_Ring) 
                gridness2=nan;
                return 
            end        
            % compute pearson correlation between rotate Acor and original Acor
            corrValues(jj) = PointCorr(org_Ring,rot_Ring);
            clear rot_Ring tmp_center_x tmp_center_y tmp_dist_from_center Y X    
        end
        % min of higher correlation at 60 and 120 degree rotation
         min_rot_60_120 = min([corrValues(3),corrValues(5)]);
        % max of lower correlations 30,90,150 rotation
        max_rot_30_90_150 = max([corrValues(2),corrValues(4),corrValues(6)]);
        % calculate gridness min(60,120)-max(30,90,150)
        gridness2 = min_rot_60_120 - max_rot_30_90_150;
    % different way to calculate gridness
    gridness1=mean(([corrValues(3),corrValues(5)]))-...
    mean([corrValues(2),corrValues(4),corrValues(6)]);
end

function out_mat=PointCorr(Rxx,RxxR)
nan_mat=Rxx .* RxxR;                                             
        notnan_inds = find(~isnan(nan_mat));  %normalized to the number of nontnan components (average)        
        n=length(notnan_inds);    
        if n < 2
            out_mat = NaN;     
        end
        sigma_x_y =sum(nan_mat(notnan_inds));
        sigma_x =      sum(Rxx(notnan_inds));
        sigma_y =      sum(RxxR(notnan_inds));
        sigma_x2 = sum(Rxx(notnan_inds).^2);
        sigma_y2 = sum(RxxR(notnan_inds).^2);
        
        out_mat = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
                                    sqrt(n*sigma_x2-sigma_x.^2) ./ ...
                                    sqrt(n*sigma_y2-sigma_y.^2);
end


 

function out_mat=Cross_Correlation(mat1,mat2)

 [ma,na] = size(mat1);
 [mb,nb] = size(mat2);
 %mc = max([ma+mb-1,ma,mb]);
 %nc = max([na+nb-1,na,nb]);
 %m=min(mc,nc);
 %out_mat = nan(m,m);
 out_mat = nan(ma+mb-1,na+nb-1);

i_size = size(mat2,1); j_size = size(mat2,2);
[work_mat,npad_i,npad_j] = pad_edges(mat1,mat2,1);

   for i = 1:size(out_mat,1)
        for j = 1:size(out_mat,2)
               
        % for each i and j, choose the correct sub-mat (size of mat 2) to
        % multiply with mat2
        
        sub_mat = work_mat(npad_i+i-floor(i_size):npad_i+i-1, ...
                                                          npad_j+j-floor(j_size):npad_j+j-1  ); 
        nan_sub_mat=sub_mat .* mat2;                                             
        notnan_inds = find(~isnan(nan_sub_mat));  %normalized to the number of nontnan components (average)
        
        n=length(notnan_inds);
        
        if n < 20
            out_mat(i,j) = NaN;
            continue;
        end
        
        sigma_x_y =sum(nan_sub_mat(notnan_inds));
        sigma_x =      sum(sub_mat(notnan_inds));
        sigma_y =      sum(mat2(notnan_inds));
        sigma_x2 = sum(sub_mat(notnan_inds).^2);
        sigma_y2 = sum(mat2(notnan_inds).^2);
        
        out_mat(i,j) = (n*sigma_x_y - sigma_x.*sigma_y) ./ ...
                                    sqrt(n*sigma_x2-sigma_x.^2) ./ ...
                                    sqrt(n*sigma_y2-sigma_y.^2);
         end % for j
    end % for i
    disp('')
end




function out_mat = nanconv2(mat,h)
    out_mat = mat;
    nan_mat = isnan(mat);
     % dilate nan_mat
    SE = strel('disk', 2);
    nan_mat =  ~imdilate(~nan_mat,SE);
    i_size = size(h,1); j_size = size(h,2);
    [work_mat,npad_i,npad_j] = pad_edges(mat,h,2);
    for i = 1:size(mat,1)
        for j = 1:size(mat,2)
            % for each i and j, choose the correct sub-mat (size of h) to multiply with h
            sub_mat = work_mat(npad_i+i-floor(i_size/2):npad_i+i+floor(i_size/2), ...
                                                              npad_j+j-floor(j_size/2):npad_j+j+floor(j_size/2)  ); % assumes h is odd in number

            notnan_inds = find(~isnan(sub_mat));  
            if ~isempty(notnan_inds)
                sum_h = sum(h(notnan_inds));   % normalize to the places without a NaN
                out_mat(i,j) = nansum(nansum(sub_mat .* h));
                out_mat(i,j) = out_mat(i,j)/sum_h;
            end
        end % for j
    end % for i
    out_mat(nan_mat) = NaN;
    disp('')
end


function [out_mat,npad_i,npad_j] = pad_edges(mat,h,l)
    npad_ij = ceil(size(h)/l);
    npad_i = npad_ij(1);
    npad_j = npad_ij(2);
    in_size = size(mat);
    out_size = in_size + [2*npad_i 2*npad_j];
    out_mat = nan(out_size);
    out_mat(npad_i+1:npad_i+in_size(1),npad_j+1:npad_j+in_size(2)) = mat;
    disp('')
end
