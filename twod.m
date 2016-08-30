clc
clear all

l = input('Enter the total length of the element: ');
b = input('Enter the total breadth of the element: ');
t = input('Enter the thickness of the element: ');

nx = input('Enter the total number of division along x axis of the element: ');
ny = input('Enter the total number of division along y axis of the element: ');

E = input('Enter the modulus of elasticity of the element: ');
mu = input('Enter the poison ratio of the element: ');

px = input('Enter the x coordinate of the point: ');
py = input('Enter the y coordinate of the point: ');

% ele = ((ny+1)*nx)+((nx+1)*ny)+(ny*nx);
% node = ((ny+1)*nx)+((nx+1)*ny);  

%% Element connectivity matrix

element(1,1)=1;
count1=0;
a = 2;
for i = 2:(2*nx)
    count1 = count1 + 1;
    if(count1==3)
        a = a +1;
        count1 = 1;
    end
    element(i,1) = a;
end

b=2;
c = nx + 2;
count1 = 0;
for i = 1:2*nx
    count1 = count1 + 1;
    if(count1==3)
       b = b + 1;
       c = c + 1;
       count1 = 1;
    end
    element(i,2) = b;
    element(i,3) = c;
end

for i = 2:2:2*nx
    element(i,2) = element(i,2) + nx + 1;
end
if(ny>1)
    for i = 2:ny
        for j = 1:8
            element(((i)*nx)+j,1) = element(((i-2)*nx)+j,1) + nx + 1;
            element(((i)*nx)+j,2) = element(((i-2)*nx)+j,2) + nx + 1;
            element(((i)*nx)+j,3) = element(((i-2)*nx)+j,3) + nx + 1;
        end
    end
end

%% nodal connectivity matrix

a = 0;
b = 0;
count1 = 0;
node = (nx+1)*(ny+1);
for i = 1:node
    nodal(i,1) = a;
    nodal(i,2) = b;
    a = a +1;
    count1 = count1 + 1;
    if(count1 ==(nx+1))
        a = 0;
        b = b + 1;
        count1 = 0 ;
    end
end
%%
% if(py - floor(py)==0)
%     p1y = floor(py) - 1;
% else
%     p1y = floor(py);
% end
% p2y = py - floor(py);
% p2x = px - floor(px);
% % if(px-floor(px)>0.5)
% %     
% %     p1x = round(px);
% %     eleno = (2*p1x)+(p1y*8);
% % else
% %     p1x = floor(px+0.49);
% %     eleno = (2*p1x) + 1+(p1y*8);
% % end
% if(p2x+p2y<=1)
%     
%     p1x = floor(px);
%     disp('huray1');
%     eleno = (2*p1x)+1+(p1y*8);
% else
% %     p1x = floor(px+0.49);
%     p1x = round(px);
%     disp('huray');
%     eleno = (2*p1x) + 2 +(p1y*8);
% end
% disp(eleno);

for i = 1:2*nx*ny
    x1 = nodal(element(i,1),1);
    y1 = nodal(element(i,1),2);
    x2 = nodal(element(i,2),1);
    y2 = nodal(element(i,2),2);
    x3 = nodal(element(i,3),1);
    y3 = nodal(element(i,3),2);
    X = [x1 x2 x3];
    Y = [y1 y2 y3];
    A = polyarea(X, Y);
    
    X1 = [px x2 x3];
    Y1 = [py y2 y3];
    A1 = polyarea(X1, Y1);
    
    X2 = [x1 px x3];
    Y2 = [y1 py y3];
    A2 = polyarea(X2, Y2);
    
    X3 = [x1 x2 px];
    Y3 = [y1 y2 py];
    A3 = polyarea(X3, Y3);
    
    if(A1 + A2 + A3 == A)
        break;
    end
end


D = E/(1-(mu^2))*[1 mu 0; mu 1 0; 0 0 0.5-(0.5*mu)];
A = 0.5;
B = 1/(2*A)*[(y2-y3) 0 (y3-y1) 0 (y1-y2) 0; 0 (x3-x2) 0 (x1-x3) 0 (x2-x1); (x3-x2) (y2-y3) (x1-x3) (y3-y1) (x2-x1) (y1-y2)];
k = t*A*B'*D*B;
disp(k);
