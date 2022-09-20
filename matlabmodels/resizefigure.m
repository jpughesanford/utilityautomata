function pos = resizefigure(scalex,scaley)

pos = get(gcf,'Position');
pos = set(gcf,'Position',[pos(1) pos(2) pos(3)*scalex pos(4)*scaley]);

end