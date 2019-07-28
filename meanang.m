function a= meanang(arads,rad)
    a=rad2deg(atan(sum(cos(arads))/sum(sin(arads))));
    if nargin==2
        a=(atan(sum(cos(arads))/sum(sin(arads))));
    end
end