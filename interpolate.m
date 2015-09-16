function y1=interpolate(xvals,yvals,x1)

% interpolates the y value corresponding to x1, assuming linearity

sizex=size(xvals);
sizey=size(yvals);

if sizex(1) == 1
    xvals=xvals';
end

if sizey(1) == 1
    yvals=yvals';
end

data = sortrows([xvals yvals]);

X = data(:,1);
Y = data(:,2);

if length(X)~=length(Y)
    error('interpolate:badinput','X and Y are of unequal length');
end

dy=diff(Y);
dx=diff(X);
minx=min(X);
maxx=max(X);

subset=find(dx>0);

for i=1:length(x1);

    search=find(X==x1(i));
    if length(search)>0
        y1(i) = mean(Y(search));
    else
        if x1 > maxx
            y1(i) = Y(length(X)) + dy(subset(length(subset)))/dx(subset(length(subset))) * (x1(i)-X(subset(length(subset))));
        else
            [search, index]= min(abs(X(subset)-x1(i)));
            y1(i) = Y(index) + dy(subset(index))/dx(subset(index)) * (x1(i)-X(subset(index)));
        end
    end
    
end