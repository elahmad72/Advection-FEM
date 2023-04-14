clc
close all
clear all
%  An Advection and diffusion equation is solved by this code
% u_x  ∂c/∂x+u_y  ∂c/∂y=D((∂^2 c)/(∂x^2 )+(∂^2 c)/(∂y^2 ))
% Quadrilateral method

%% enter the number of element in each side ( square shape is considered)
N1=input('enter number of elements in each row');

N=power(N1,2);
%% the number of nodes is calculated
node_num=power(N1+1,2);

%% elem matrix is the definition of each element according to it's nodes
elem=ones(N,4);
dx=1/N1;
for i=1:N
    for j=1:4
        if j==1
            if mod(i,N1)==0
                elem(i,j)=floor(i/(N1))+i-1;
            else
                elem(i,j)=floor(i/(N1))+i;
            end
        elseif j==2
            elem(i,j)=elem(i,j-1)+1;
        elseif j==3
            elem(i,j)=elem(i,j-1)+N1;
        else
            elem(i,j)=elem(i,j-1)+1;
        end
    end
end
%% node matrix is the location of each element according to the xy-coordination
node=zeros(node_num,2);

for i=1:node_num
    for j=1:2
        if j==1
            if mod(i,N1+1)==0
                node(i,j)=1;
            else
                node(i,j)=(mod(i,(N1+1))-1)*1/N1;
            end
        elseif j==2
            if mod(i,N1+1)==0
                node(i,j)=(floor(i/(N1+1))-1)*1/N1;
            else
                node(i,j)=floor(i/(N1+1))*1/N1;
            end
            
        end
    end
end
%%  bdFlag matrix is the state of boundary condition each element;
%  0 for non-boundary sides;
%  1 for the first type, i.e., Dirichlet boundary;
%  2 for the second type, i.e., Neumann boundary;
%  3 for the third type, i.e., Robin boundary
bdFlag=zeros(N1+1,4);

for i=1:N
    for j=1:4
        if j==1
            if i/(N1)<=1
                bdFlag(i,j)=1;
            end
        elseif j==2
            if mod(i,N1)==0
                bdFlag(i,j)=1;
            end
        elseif j==3
            if i/(N1)>N1-1
                bdFlag(i,j)=1;
            end
        else
            if mod(i,N1)==1
                bdFlag(i,j)=1;
            end
            
        end
    end
end

% % %  Calculation
% constants
D= 8*10^-5;   %mass transfer coeficient
u_x=10;     %x-vector of velocity
u_y=10;      %y-vector of velocity
c0= .05;       % saturated concentration in Dirichlet boundary

% %  matrix of global weight

% % % % % % Matrix of local weight

b=dx/2;
c=dx/2;
w44=zeros(4,4);

for i=1:4
    w44(i,i)=-((u_x-D)*b^2+(u_y-D)*c^2)/(3*b*c);
end
w44(1,2)=-((u_x-D)*b^2-2*(u_y-D)*c^2)/(6*b*c);
w44(1,3)=((u_x-D)*b^2+(u_y-D)*c^2)/(6*b*c);
w44(1,4)= -((u_x-D)*c^2-2*(u_y-D)*b^2)/(6*b*c);
w44(3,4)=w44(1,2);
w44(2,3)=w44(1,4);
w44(2,4)=w44(1,3);
w44=w44+triu(w44,1)';

% % % % % % Matrix of global wieght
wglobe=zeros(node_num);
for i=1:node_num-3
    for j=1:4
        for k=1:4
            wglobe(i+j-1,i+k-1)=wglobe(i+j-1,i+k-1)+w44(j,k);
        end
    end
end

%% Matrix of Concentration
C=ones(node_num,1);
[a,b]=find(node(:,1)==0);
[aa,bb]=find(node(:,2)==0);
[aaa,bbb]=find(node(:,2)==1);
[aaaa,bbbb]=find(node(:,1)==1);
for i=a
    C(i,1)=c0;
end

for i=aa
    C(i,1)=.01;
end
for i=aaa
    C(i,1)=.01500;
end
for i=aaaa
    C(i,1)=c0/2;
end
node1=node;
for i=1: length(a)
    wglobe(:,a(i))=0;
    wglobe(a(i),:)=0;
    wglobe(a(i),a(i))=1;
    wglobe(:,aa(i))=0;
    wglobe(aa(i),:)=0;
    wglobe(aa(i),aa(i))=1;
    wglobe(:,aaa(i))=0;
    wglobe(aaa(i),:)=0;
    wglobe(aaa(i),aaa(i))=1;
    wglobe(:,aaaa(i))=0;
    wglobe(aaaa(i),:)=0;
    wglobe(aaaa(i),aaaa(i))=1;
end

%% slover
CC=linsolve(wglobe,C);


%% Calculation for concentration of elements
AAA = [a; aa;aaa;aaaa];
AAA=unique(AAA);
A(:,1)=1:node_num;
AA=setdiff(A,AAA);

C_elem=zeros(N,1);
error=2;
Tol=0.0001;
counter=0;
C1=C;
while sum(abs(error))>Tol
    CC=linsolve(wglobe,C1);
    
    for i=1:N
        S(i,1)=1/(4*b*c)*(b-node(elem(i,1),1))*(c-node(elem(i,1),2));
        S(i,2)=1/(4*b*c)*(b+node(elem(i,2),1))*(c-node(elem(i,2),2));
        S(i,3)=1/(4*b*c)*(b+node(elem(i,3),1))*(c+node(elem(i,3),2));
        S(i,4)=1/(4*b*c)*(b-node(elem(i,4),1))*(c+node(elem(i,4),2));
        C_elem(i,1)=S(i,:)*[CC(elem(i,1)); CC(elem(i,2)); CC(elem(i,3)); CC(elem(i,4))];
    end
    
    
    for i=1:length(AA)
        
        [tt,b]=find(elem(:,:)==AA(i));
        for j=1:length(tt)
            C1(AA(i),1)=CC(AA(i),1)+1/length(tt)*C_elem(j,1);
            
        end
    end
    error=C1-CC;
    
    counter=counter+1;
    
end

%% Printing figures
%%3D figure
for i=1:N1+1
    C1(1:N1+1,i)=C1((i-1)*(N1+1)+1:(i-1)*(N1+1)+N1+1,1);
    
    
end
C1(N1+2:node_num,:)=[];
figure(1)
hold on
[X,Y] = meshgrid(0:1/(N1):1);
surf(Y,X,C1)

hold off
