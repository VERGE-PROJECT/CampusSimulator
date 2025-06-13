function [r,inters_point] = intersect_lines(xt,yt,xr,yr,x1,y1,x2,y2)
%Returns r=1 if the line (xt,yt)->(xr,yr) intersects the line
%(x1,y1)->(x2,y2)
%It also returns the intersection point inters_point(x,y)

if(x2==x1)
    y=yt+(yr-yt)*(x1-xt)/(xr-xt);
    x=x2;
    if (y>=min(yt,yr))&&(y<=max(yt,yr))&&(y>=min(y1,y2))&&(y<=max(y1,y2))
        r=1;
    else
        r=0;
    end
elseif (xt==xr)
    y=y2+(y2-y1)*(xr-x2)/(x2-x1);
    x=xt;
    if (y>=min(yt,yr))&&(y<=max(yt,yr))&&(y>=min(y1,y2))&&(y<=max(y1,y2))
        r=1;
    else
        r=0;
    end
else
    x=((yt-y2)*(x2-x1)*(xr-xt)+(y2-y1)*x2*(xr-xt)-(yr-yt)*xt*(x2-x1))/((y2-y1)*(xr-xt)-(x2-x1)*(yr-yt));
    y=y2+(x-x2)*(y2-y1)/(x2-x1);
    if (x>=min(xt,xr))&&(x<=max(xt,xr))&&(x>=min(x1,x2))&&(x<=max(x1,x2))
        r=1;
    else
        r=0;
    end
end
inters_point=[x,y];




