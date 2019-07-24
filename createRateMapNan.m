function [rm, max_r] = createRateMapNan(c, nbins)
px=    c.px;
py=    c.py;
pt=    c.pt;
sx=    c.sx;
sy=    c.sy;
st=    c.st;
if ~isempty(st)
    px = toCol(px); py = toCol(py); %pt = toCol(pt);
    sx = toCol(sx); sy = toCol(sy); %st = toCol(st);
    %si = discretize(st, [-Inf; mean([pt(2:end) pt(1:end-1)],2); +Inf]);
    mx = max(px);
    my = max(py);
    pxi = discretize(px, 0:mx/nbins:mx);
    pyi = discretize(py, 0:my/nbins:my);
    rmt = accumarray([pyi, pxi], 1, [nbins nbins]);
    %rmt(rmt<1) = 1e6; %NO NANS??
    sxi = discretize(sx, 0:mx/nbins:mx); %px(si)
    syi = discretize(sy, 0:my/nbins:my); %py(si)
    rms = accumarray([syi,sxi], 1, [nbins nbins]);
    %!!!!
    %rm = rms./(rmt*.02);
    %max_r = max(rm(:));
    %add in extra timestep ONLY when spikes occured faster than timestep
    %t = rms-rmt; t(t<0)=0; rmt = rmt + t; 
    %if sum(t)~=0
    %    disp('$$ STILL ADDING TIME TO BINS');
    %end
    %max_r = (max(t(:))+1)/.02;
    rmt(rmt==0&rms~=0)=1; %add to bin where 
    rmt=rmt*0.02;
    rm = rms ./ rmt;
    max_r = max(rm);
    %rate_mat(isnan(rate_mat)) = 0; %DO THIS?
else
    rm = nan(nbins);
    max_r = 0;
end


end
