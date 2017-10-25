function [p,e] = surfaceIsocontour(f,v,data,isoval);
if nargin < 3
    error('Must have at least 3 input arguments');
end
if nargin < 4
    isoval = max(data(:)) + min(data(:))/2.0;
end

% clear all;
% close all;
% addpath /cis/home/dtward/Functions/byu
% [f,v] = readbyu_local('/cis/home/dtward/Documents/Mostofsky/atlas_2693/atlas_2693_01.byu');
% data =  (v(:,1)/max(abs(v(:,1))))  + (v(:,2)/max(abs(v(:,2)))).^2 + (v(:,3)/max(abs(v(:,3)))).^3;
% data = sin(2*pi*v(:,1)/10)
% isoval = (max(data(:))+min(data(:)))/2.0;
% h = patch('faces',f,'vertices',v,'facevertexcdata',data,'facecolor','interp','facelighting','gouraud');
% axis image
% light

% now we loop through each triangle
p = zeros(size(f,1)*2,3); % maximum size
count = -1;
for floop = 1 : size(f,1)
    % check the value of the function at the three vertices
    
    v1 = v(f(floop,1),:);
    v2 = v(f(floop,2),:);
    v3 = v(f(floop,3),:);
    
    d1 = data(f(floop,1));
    d2 = data(f(floop,2));
    d3 = data(f(floop,3));
    
    b1 = d1>isoval;
    b2 = d2>isoval;
    b3 = d3>isoval;
    b = [b1,b2,b3];
%     disp(b)
    
    % now we go through each case
    if prod(b == [0, 0, 0]) % nobody is above the threshold
%         disp('000')
        continue
    elseif prod(b == [0,0,1]) % point 3 is above the threshold
        count = count + 2;
        % we draw a line from edge 31, to edge 32
        t = abs(d3-isoval)/abs(d3-d1);
        p1 = v3*(1-t) + v1*t;
        t = abs(d3-isoval)/abs(d3-d2);
        p2 = v3*(1-t) + v2*t;
        p(count + (0:1),:) = [p1; p2];
    elseif prod(b == [0,1,0]) % point 2 is above the threshold
        count = count + 2;
        % we draw a line from edge 21 to 23
        t = abs(d2-isoval)/abs(d2-d1);
        p1 = v2*(1-t) + v1*t;
        t = abs(d2-isoval)/abs(d2-d3);
        p2 = v2*(1-t) + v3*t;
        p(count + (0:1),:) = [p1; p2];
    elseif prod(b == [0,1,1])
        count = count + 2;
        % we draw a line from edge 12 to edge 13
        t = abs(d1-isoval)/abs(d1-d2);
        p1 = v1*(1-t) + v2*t;
        t = abs(d1-isoval)/abs(d1-d3);
        p2 = v1*(1-t) + v3*t;
        p(count + (0:1),:) = [p1; p2];        
    elseif prod(b == [1,0,0])
        count = count + 2;
        % we draw a line from edge 12 to 13
        t = abs(d1-isoval)/abs(d1-d2);
        p1 = v1*(1-t) + v2*t;
        t = abs(d1-isoval)/abs(d1-d3);
        p2 = v1*(1-t) + v3*t;
        p(count + (0:1),:) = [p1; p2];    
    elseif prod(b == [1,0,1])
        count = count + 2;
        % we draw a line from 21 to 23
        t = abs(d2-isoval)/abs(d2-d1);
        p1 = v2*(1-t) + v1*t;
        t = abs(d2-isoval)/abs(d2-d3);
        p2 = v2*(1-t) + v3*t;
        p(count + (0:1),:) = [p1; p2];    
    elseif prod(b == [1,1,0]) % point 3 is below the threshold
        count = count + 2;
        % we draw a line from edge 31, to edge 32
        t = abs(d3-isoval)/abs(d3-d1);
        p1 = v3*(1-t) + v1*t;
        t = abs(d3-isoval)/abs(d3-d2);
        p2 = v3*(1-t) + v2*t;
        p(count + (0:1),:) = [p1; p2];
    elseif prod(b == [1,1,1]) % nobody is below, don't draw anything
%         disp('111');
        continue
    end
    
end
p = p(1:(count+1),:);

% figure;
e = 1:size(p,1);
e = reshape(e,2,[]);
e = e';
if nargout == 0
    patch('vertices',p,'faces',e,'edgecolor','k','linewidth',2,'facecolor','none')
end
