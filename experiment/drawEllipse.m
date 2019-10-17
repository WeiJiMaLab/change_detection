% Output an ellipse with axes length d1 and d2, rotation rot. fgcol is the
% foreground color, bgcol is the background color. im is a double matrix.

function im = drawEllipse(d1,d2,rot,fgcol,bgcol)

rot = -rot-90;  

d1 = round(d1);
d2 = round(d2);

rot = -rot/180*pi;

% make sure that d1 is the minor axis
if (d1>d2)
    d3=d1;
    d1=d2;
    d2=d3;
end

% draw ellipse
im = ones(2*d2,2*d2)*bgcol;
minX = -d2;
maxX = minX + 2*d2 - 1; 
[X Y] = meshgrid(minX:maxX,minX:maxX);
X_new = X * cos(rot) - Y * sin(rot);
Y = X * sin(rot) + Y * cos(rot);
X = X_new;
idx = (X.^2/(d1/2)^2 + Y.^2/(d2/2)^2)<1;
im(idx) = fgcol;

% crop
while im(:,1)==bgcol
    im = im(:,2:end);
end
while im(1,:)==bgcol
    im = im(2:end,:);
end
while im(end,:)==bgcol
    im = im(1:end-1,:);
end
while im(:,end)==bgcol
    im = im(:,1:end-1);
end


