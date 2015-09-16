function y=SmoothedLinearBaselineRemove(fields,dat,xpoints)

y.fields=fields;
y.DatSmoothed=smooth(dat,21);
y.fit=fit(fields,dat,'linearinterp');
linearbaseline=fit([xpoints(1) ; xpoints(2)],[y.fit(xpoints(1)) ; y.fit(xpoints(2)) ],'linearinterp');
y.subtracted=y.DatSmoothed-linearbaseline(fields);
y.subtractedfit=fit(fields,y.subtracted,'linearinterp');
y.AreaUnder=integrate(y.subtractedfit,xpoints(2),xpoints(1));
y.PeakMax=max(y.fit(xpoints(1):(xpoints(2)-xpoints(1))/100:xpoints(2)));
y.AreaUnderToPeakMax=y.AreaUnder/y.PeakMax;