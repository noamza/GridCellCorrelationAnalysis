function groups = find_simultaneously_recorded_cells(cells)
    groups = [];
    for i=1:length(cells)
        a = cells{i};
        if ~isempty(a.middle) %skip cells with non middle
            key = sprintf('g_%s_%s', a.id, a.date);
            if isfield(groups, key)
                t = groups.(key);
                t = [t a];
                groups.(key) = t;
            else
                groups = setfield(groups, key, a);
            end
        end
        %    fprintf('%d: bad cell %f rmthr %f', i, a.max_r, rate_map_thresh);
    end
    
    %validate
    keys = fieldnames(groups);
    for i = 1:length(keys)
        k = keys{i};
        group = groups.(k);
        t = group(1);
        t = length(t.before.px);
        for j = 1:length(group) %CHECK TIMES !!
            tt = group(j); %cells{group(j)};
            tt = length(tt.before.px);
            if t ~= tt
                fprintf('%f %f whats up with this group b.length %s?\n',t, tt, k);
            end
        end
    end
    
    fprintf('grouped\n');
end