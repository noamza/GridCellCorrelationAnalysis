function a = patchTrajectoryLinear(pt,px,py,dt,dtthresh)
    pt = round(pt,3); dt = round(dt,3); dtthresh = round(dtthresh,3);
    if isnan(dt)
            a.x = 0; 
            a.y = 0; 
            a.t = 1;
        return
    end
    t = zeros(ceil(pt(end)/dt),1); x = t; y = t;
    t(1) = pt(1); x(1) = px(1); y(1) = py(1);
    ii = 1;
    for i = 2:length(pt)
       if pt(i) - pt(i-1) > dtthresh %%pt(i) - pt(i-1)
            tt = pt(i-1):dt:pt(i);
            tt = tt(2:end-1);
            step = dt/(pt(i)- pt(i-1)); %number < 1
            assert(abs(step)<1);
            r = (1:length(tt))*step;
            tx = px(i-1) + r*(px(i)- px(i-1));
            ty = py(i-1) + r*(py(i)- py(i-1));
            %gap not large enough to add timestep
            if(isempty(tt))
                tx = px(i);
                ty = py(i);
                tt = pt(i);
            end
       else
            tx = px(i);
            ty = py(i);
            tt = pt(i);
       end
       ii = ii+length(tt);
       t(ii-length(tt)+1:ii) = tt;%t = [t tt];
       x(ii-length(tt)+1:ii) = tx;%x = [x tx];
       y(ii-length(tt)+1:ii) = ty;%y = [y ty];
    end
    a.x = x(1:ii); 
    a.y = y(1:ii); 
    a.t = t(1:ii);
    assert(min(diff(a.t))>dt/2);
    assert(max(diff(a.t))<2*dtthresh,'max(diff(t))<dtthresh');
    assert(sum(a.t==0)<=1,'sum(isnan(a.x))==0 trajectory');
    
%     is = 1:length(pt);
%     t = diff(pt);
%     is = is(t>dtthresh);
%     patches = cell(1,length(is));
%     for i = 1:length(is)
%         %i = 441
%         t = pt(is(i)):dt:pt(is(i)+1);
%         patches(i).t = t(2:end-1);
%         step = dt/(pt(is(i)+1)- pt(is(i))) %number < 1
%         assert(abs(step)<1);
%         r = (1:length(t))*step;
%         patches(i).x = px(is(i)) + r*(px(is(i)+1)- px(is(i)));
%         patches(i).y = py(is(i)) + r*(py(is(i)+1)- py(is(i)));
%         [ c.px(is(i)) x c.px(is(i)+1)]
%         [ c.py(is(i)) y c.py(is(i)+1)]
%         [ c.pt(is(i)) t c.pt(is(i)+1)]
%         diff([ c.px(is(i)) x c.px(is(i)+1)])
%         diff([ c.py(is(i)) y c.py(is(i)+1)])
%         diff([ c.pt(is(i)) t c.pt(is(i)+1)])
    
end

