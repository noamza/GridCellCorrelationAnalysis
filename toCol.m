function x = toCol(x)
    if ~iscolumn(x)
        x = x';
    end
end

function x = toRow(x)
    if iscolumn(x)
        x = x';
    end
end
