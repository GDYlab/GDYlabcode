function tbins = get_tbins(ts, te, tau, percent_overlap)
    tbins = ts:tau*(1-percent_overlap):te;
    if numel(tbins)==1
        tbins = [tbins, tbins+tau];
    end
    tbins = [tbins', tbins'+tau];
    
    if tbins(end,2) > te
        tbins = tbins(1:end-1,:);
    end
end