% Central Catadioptric Camera Calibration Using Line Image and Mirror Contour
clc
clear
close all

%%
%Basic parameters%
K=[5000 0.1  2000; %Intrinsic parameters, K,3*3，5 degrees of freedom%
    0   4000 1000;
    0   0    1];
l=0.2; %Mirror parameter, (0,1)%

%%
%Step 1, get image%

%Mirror Contour%
% Q0=[1 0 0;
%     0 1 0;
%     0 0 -1/(l^2)]; %1/l为tan值，为理想镜面，视角为180°

d0=0; %镜面平面与单位视球赤道的垂直距离
Q0=[1 0 0;
    0 1 0;
    0 0 -(1-d0^2)/(l+d0)^2]; %sqrt(1-d0^2)/(l+d0)为tan值，为非理想镜面，视角小于180°
    
CQ0=inv(K')*Q0*inv(K);
CQ0=CQ0/CQ0(3,3);
cQ0=[CQ0(1,1) 2*CQ0(1,2) CQ0(2,2) 2*CQ0(1,3) 2*CQ0(2,3) CQ0(3,3)]';


%Line image%
x1=[-4;5;5];
x2=[3;6;5];

n1=x1/norm(x1);
d1=0;
n2=x2/norm(x2);
d2=0;


Q1=[(l^2-1)*n1(1)^2+(d1+l*n1(3))^2  (l^2-1)*n1(1)*n1(2)                -(l*d1+n1(3))*n1(1);
    (l^2-1)*n1(1)*n1(2)               (l^2-1)*(n1(2))^2+(d1+l*n1(3))^2   -(l*d1+n1(3))*n1(2);
    -(l*d1+n1(3))*n1(1)                -(l*d1+n1(3))*n1(2)                 d1^2-n1(3)^2      ];
% Q1=Q1/Q1(3,3);

Q2=[(l^2-1)*n2(1)^2+(d2+l*n2(3))^2  (l^2-1)*n2(1)*n2(2)                -(l*d2+n2(3))*n2(1);
    (l^2-1)*n2(1)*n2(2)               (l^2-1)*(n2(2))^2+(d2+l*n2(3))^2   -(l*d2+n2(3))*n2(2);
    -(l*d2+n2(3))*n2(1)                -(l*d2+n2(3))*n2(2)                 d2^2-n2(3)^2      ];


C1=inv(K')*Q1*inv(K);
C1=C1/C1(3,3);
c1=[C1(1,1) 2*C1(1,2) C1(2,2) 2*C1(1,3) 2*C1(2,3) C1(3,3)]';


C2=inv(K')*Q2*inv(K);
C2=C2/C2(3,3);
c2=[C2(1,1) 2*C2(1,2) C2(2,2) 2*C2(1,3) 2*C2(2,3) C2(3,3)]';


%%
%Step 2, find vertices%

[V1,e1]=eig(C1,CQ0);
e1=diag(e1);
if isreal(e1) 
    [~,k1]=max(abs(median(sign(e1))-sign(e1))); 
else 
    for i=1:3
        if isreal(e1(i))
            k1=i;
        end
    end
end
v1=V1(:,k1); %get vertex v1
l01=C1*v1; %get the line passing through the principal point


[V2,e2]=eig(C2,CQ0);
e2=diag(e2);
if isreal(e2) 
    [~,k2]=max(abs(median(sign(e2))-sign(e2)));
else 
    for i=1:3
        if isreal(e2(i))
            k2=i;
        end
    end
end
v2=V2(:,k2); %get vertex v2
l02=C2*v2; %get the line passing through the principal point


%%
%Step 3, obtain the simplified ICPs%

Op=cross(l01,l02); 
Op=Op/Op(3);

Tp=[1 0 0;
    0 1 0;
    Op']';

C1p=Tp'*C1*Tp; 
C2p=Tp'*C2*Tp; 

l1=C1p*[0 0 1]';
l2=C2p*[0 0 1]';

%
%Function, find the intersections between line L[A1;B1;D1] and conic C=[a1 b1 c1 d1 e1 1]'%

function [X,Y]=crossLandC(L,C)
A1=L(1);
B1=L(2);
D1=L(3);

C=C/C(3,3);
a1=C(1,1);
b1=2*C(1,2);
c1=C(2,2);
d1=2*C(1,3);
e1=2*C(2,3);

X=zeros(2,1);
Y=zeros(2,1);

X(1)=-(D1 + (B1*(A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1^2*conj(e1) + A1*B1*conj(d1) + A1*D1*conj(b1) - 2*B1*D1*conj(a1)))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2)))/A1;

X(2)=-(D1 - (B1*(A1^2*conj(e1) + A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1*B1*conj(d1) - A1*D1*conj(b1) + 2*B1*D1*conj(a1)))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2)))/A1;

Y(1)=(A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1^2*conj(e1) + A1*B1*conj(d1) + A1*D1*conj(b1) - 2*B1*D1*conj(a1))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2));

Y(2)=-(A1^2*conj(e1) + A1*(A1^2*conj(e1)^2 + B1^2*conj(d1)^2 + D1^2*conj(b1)^2 - 4*B1^2*conj(a1) - 4*A1^2*conj(c1) - 4*D1^2*conj(a1)*conj(c1) + 4*A1*B1*conj(b1) - 2*A1*B1*conj(d1)*conj(e1) - 2*A1*D1*conj(b1)*conj(e1) + 4*A1*D1*conj(c1)*conj(d1) + 4*B1*D1*conj(a1)*conj(e1) - 2*B1*D1*conj(b1)*conj(d1))^(1/2) - A1*B1*conj(d1) - A1*D1*conj(b1) + 2*B1*D1*conj(a1))/(2*(conj(c1)*A1^2 - conj(b1)*A1*B1 + conj(a1)*B1^2));

end

%

[X1,Y1]=crossLandC(l1,C1p);
mi1=[X1(1),Y1(1),1]';
mj1=[X1(2),Y1(2),1]';

[X2,Y2]=crossLandC(l2,C2p);
mi2=[X2(1),Y2(1),1]';
mj2=[X2(2),Y2(2),1]';

%%
%Step 4, obtain the intrinsic and mirror parameters%

A=[real(mi1(1)^2),real(mi1(1)*mi1(2)),real(mi1(2)^2),1;
   imag(mi1(1)^2),imag(mi1(1)*mi1(2)),imag(mi1(2)^2),0;
   real(mi2(1)^2),real(mi2(1)*mi2(2)),real(mi2(2)^2),1;
   imag(mi2(1)^2),imag(mi2(1)*mi2(2)),imag(mi2(2)^2),0];

[~,~,V]=svd(A); 
c=V(:,end);
c=c/c(4);
C=[c(1)   c(2)/2 0;
   c(2)/2 c(3)   0;
   0      0      c(4)];
K1=inv(chol(C)); 
K1=K1/K1(3,3) % get K

Qs1=KKK'*C1*KKK; 
Qs1=Qs1/Qs1(3,3);
Qs2=KKK'*C2*KKK; 
Qs2=Qs2/Qs2(3,3);

l1A=sqrt(1-(Qs1(3,3)*Qs1(2,1))/(Qs1(3,2)*Qs1(3,1))); 
l2A=sqrt(1-(Qs2(3,3)*Qs2(2,1))/(Qs2(3,2)*Qs2(3,1))); 

lAA=(l1A+l2A)/2  % get l


