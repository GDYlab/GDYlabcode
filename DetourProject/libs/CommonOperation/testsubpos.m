function testsubpos(pos)

if iscell(pos)
    for isub = 1:length(pos)
        ax = subplot('position',pos{isub});
        box(ax,'on')
        xticks(ax,[])
        yticks(ax,[])
    end
else
    ax = subplot('position',pos);
    box(ax,'on')
    xticks(ax,[])
    yticks(ax,[])
end
end

