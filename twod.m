clc
clear all

l = input('Enter the total length of the element: ');
%b = input('Enter the total breadth of the element: ');
t = input('Enter the thickness of the element: ');

nx = input('Enter the total number of division along x axis of the element: ');
ny = input('Enter the total number of division along y axis of the element: ');

E = input('Enter the modulus of elasticity of the element: ');
mu = input('Enter the poison ratio of the element: ');

px = input('Enter the x coordinate of the point: ');
py = input('Enter the y coordinate of the point: ');

%% Boundary condition
nf = input('Enter the number of node fixed: ');
q = ones(2*(nx+1)*(ny+1),1);
for i = 1:nf
    nfn = input('Enter the node number which is fixed: ');
    q((2*nfn -1),1) = 0;
    q((2*nfn),1) = 0; 
end

%% Loading
nl = input('Enter the number of load acting: ');
load = zeros(2*(nx+1)*(ny+1),1);
for i = 1:nl
    nnla = input('Enter the node where load is acting: ');
    load((2*nnla)-1,1) = input('Enter the amount of horizontal force acting at that node ');
    load(2*nnla,1) = input('Enter the amount of vertical force acting at that node ');
end

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
nodal= [];
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

%% global stiffness matrix
totnode = (nx+1)*(ny+1);
temp = k;
K = zeros(totnode);

for i = 1:size(element,1)
    m(1) = (2*element(i,1))-1;
    m(2) = (2*element(i,1));
    m(3) = (2*element(i,2))-1;
    m(4) = (2*element(i,2));
    m(5) = (2*element(i,3))-1;
    m(6) = (2*element(i,3));
    for r = 1:6
        for c = 1:6
            K(m(r),m(c)) = K(m(r)+m(c)) + temp(r,c);
        end
    end
end       

%% Displacement of each node
count1 = 1;
%K1 = K ;
%load1 = load ;
cnt = 1;
for i = 1:2*totnode
    if(q(i)==0)
        inde(count1,1) = i; 
        %K1(i,:) = 0;
        %load1(i,1) = 0;
    else
        K1(cnt,:) = K(i,:);
       % K1(:,cnt) = K(:,i);
        load1(cnt,:) = load(i,:);
        cnt = cnt +1;
    end
    count1 =count1 + 1;
    %K1(i,:) = K(i,:);
    %load1(i,1) = load(i,1);
end
cnt = 1;
for i = 1:2*totnode
    if(q(i)==0)
       
    else
        K2(:,cnt) = K1(:,i);
       % K1(:,cnt) = K(:,i);
        cnt = cnt +1;
    end

end
q1 = q;
q1( ~any(q1,2), : ) = [];  % to remove element of a matix with all zero coloumn;
inde( ~any(inde,2), : ) = [];
%load1( ~any(load1,2), : ) = [];

q1 = (inv(K2))*load1;
% cnt=1;
% q2 = [];
% cnt1 = 1;
% for i = 1: 2*totnode
%     if(size(inde)>=cnt) 
%         if(inde(cnt)==i)
%             q2(i)= 0;
%             disp(['cnt is ', num2str(cnt)]);
%             cnt = cnt+ 1;
%             disp(['i is ', num2str(i)]);
%         end
%     elseif(inde(cnt)~=i)
%         disp(['i is ', num2str(i)]);
%         q2(i) = q1(cnt1);
%         disp(['cnt1 is ', num2str(cnt1)]);
%         cnt1 = cnt1 + 1;
%         disp(q1(cnt1));
%     end
%     
% end
a = inde;
cnta = 1;
cntq1 = 1;
q2 =[];
for i = 1 : size(q1,1) + size(a,1)
    if (cnta <= size(a,1) && i == a(cnta))
        q2 = [q2 ; 0];
        cnta = cnta +1;
    else
        q2 = [q2 ; q1(cntq1)];
        cntq1 = cntq1 + 1;
        
    end
end
disp(q2);

%% Strain

for i = 1:2*nx*ny
    X1 = nodal(element(i,1),1);
    Y1 = nodal(element(i,1),2);
    X2 = nodal(element(i,2),1);
    Y2 = nodal(element(i,2),2);
    X3 = nodal(element(i,3),1);
    Y3 = nodal(element(i,3),2);
    st1 = (Y2-Y3)*q2((2*element(i,1))-1)+(Y3-Y1)*q2((2*element(i,2))-1)+(Y1-Y2)*q2((2*element(i,3))-1);
    st2 = (X3-X2)*q2(2*element(i,1))+(X1-X3)*q2(2*element(i,2))+(X2-X1)*q2(2*element(i,3));
    st3 = (X3-X2)*q2((2*element(i,1))-1)+(Y2-Y3)*q2(2*element(i,1))+(X1-X3)*q2((2*element(i,2))-1)+(Y3-Y1)*q2(2*element(i,2))+(X2-X1)*q2((2*element(i,3))-1)+(Y1-Y2)*q2(2*element(i,3));

    %strain(i) = 1/(2*A)[(Y2-Y3)*q2((2*element(i,1))-1)+(Y3-Y1)*q2((2*element(i,2))-1)+(Y1-Y2)*q2((2*element(i,3))-1); (X3-X2)*q2(2*element(i,1))+(X1-X3)*q2(2*element(i,2))+(X2-X1)*q2(2*element(i,3)); (X3-X2)*q2((2*element(i,1))-1)+(Y2-Y3)*q2(2*element(i,1))+(X1-X3)*q2((2*element(i,2))-1)+(Y3-Y1)*q2(2*element(i,2))+(X2-X1)*q2((2*element(i,3))-1)+(Y1-Y2)*q2(2*element(i,3))];
    strain(:,i) = [st1; st2;st3];
    %strain(2,i) = st2;
    %strain(3,i) = st3;
    
    stress(:,i) = D*strain(:,i);
end
    