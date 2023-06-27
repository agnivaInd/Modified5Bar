clear all
close all
clc

l1=50;l2=40;l3=50;l4=50;l5=40;
phi=pi/6;

xf = [];
yf = [];
step = 0.2;

j=15;
while j<70
    j
    xf = [xf;j];
    yf = [yf;80];
    j = j + step;
end

j=80;
while j>60
    j
    xf = [xf;70];
    yf = [yf;j];
    j = j - step;
end

j=60;
while j>15
    j
    xf = [xf;j];
    yf = [yf;60];
    j = j - step;
end

j=60;
while j<80
    j
    xf = [xf;15];
    yf = [yf;j];
    j = j + step;
end

theta1 = [];
for j = 1:length(xf)
    j
    A = 2*xf(j)*l1;
    B = 2*yf(j)*l1;
    C = xf(j)^2 + yf(j)^2 + l1^2 - 4*l2^2*(cos(phi))^2;
    theta1 = [theta1;2*atan((B+sqrt(B^2+A^2-C^2))/(A+C))];
end

theta2=[];
for j = 1:length(xf)
    j
    k = (xf(j)-l1*cos(theta1(j)))/(2*l2*cos(phi));
    theta2 = [theta2;2*atan(sqrt((-k+1)/(k+1)))-phi];
end

theta4=[];
for j = 1:length(xf)
    j
    B1 = -2*l1*l4*sin(theta1(j))-2*l2*l4*sin(theta2(j));
    A1 = 2*l4*l5-2*l1*l4*cos(theta1(j))-2*l2*l4*cos(theta2(j));
    C1 = l1^2+l2^2+l4^2+l5^2-l3^2 + 2*l1*l2*cos(theta1(j)-theta2(j)) - 2*l1*l5*cos(theta1(j)) - 2*l2*l5*cos(theta2(j));
    disp((-B1-sqrt(B1^2+A1^2-C1^2))/(C1-A1));
    theta4 = [theta4;2*atan((-B1-sqrt(B1^2+A1^2-C1^2))/(C1-A1))];
end

theta3=[];
for j = 1:length(xf)
    j
    K1=(l1*sin(theta1(j))+l2*sin(theta2(j))-l4*sin(theta4(j)))/l3;
    theta3 = [theta3;2*atan((1+sqrt(1-K1^2))/K1)];
end

pA = [0 0];
pB = [l1*cos(theta1) l1*sin(theta1)];
pE = [l5 0];
pD = [l5+l4*cos(theta4) l4*sin(theta4)];
pC = [l5+l4*cos(theta4)+l3*cos(theta3) l4*sin(theta4)+l3*sin(theta3)];


grid ON
hold on

for i = 1:length(xf)

    F = [xf(i) yf(i)];
    p1=viscircles([pA(1) pA(2)],0.6);
    p2=viscircles([pB(i,1) pB(i,2)],0.6);
    p3=viscircles([pC(i,1) pC(i,2)],0.6);
    p4=viscircles([pD(i,1) pD(i,2)],0.6);
    p5=viscircles([pE(1) pE(2)],0.6);

    xs = [pA(1) pB(i,1) pC(i,1) pD(i,1) pE(1) pA(1) pB(i,1) xf(i) pC(i,1)];
    ys = [pA(2) pB(i,2) pC(i,2) pD(i,2) pE(2) pA(2) pB(i,2) yf(i) pC(i,2)];

    Fx(i) = F(1);
    Fy(i) = F(2);
    Plot=plot(xs,ys,'k','Linewidth',3);
   
    Plot2=plot(Fx,Fy,'-r','Linewidth',1);
    set(gca,'Xlim',[-50,120],'Ylim',[-50,120]);
   
    pause(0.01);
    if i<length(xf)
        delete(p1);
        delete(p2);
        delete(p3);
        delete(p4);
        delete(p5);
        delete(Plot);
        delete(Plot2);
    end
end
