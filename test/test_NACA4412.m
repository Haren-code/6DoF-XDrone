clear, clc
close all
%%% Point Style

A = [1.0000     0.0013
0.9500     0.0147
0.9000     0.0271
0.8000     0.0489
0.7000     0.0669
0.6000     0.0814
0.5000     0.0919
0.4000     0.0980
0.3000     0.0976
0.2500     0.0941
0.2000     0.0880
0.1500     0.0789
0.1000     0.0659
0.0750     0.0576
0.0500     0.0473
0.0250     0.0339
0.0125     0.0244
0.0000     0.0000
0.0125     -0.0143
0.0250     -0.0195
0.0500     -0.0249
0.0750     -0.0274
0.1000     -0.0286
0.1500     -0.0288
0.2000     -0.0274
0.2500     -0.0250
0.3000     -0.0226
0.4000     -0.0180
0.5000     -0.0140
0.6000     -0.0100
0.7000     -0.0065
0.8000     -0.0039
0.9000     -0.0022
0.9500     -0.0016
1.0000     -0.0013
1.0000     0.0013];


%%%% Equation Style
N_span  = 20;
N_air   = 120;
c_blade = 1;   % chord length (Z direction)
m = 0.04; p = 0.4; t = 0.12;
x = linspace(0,1,N_air);

yt = 5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2+0.2843*x.^3-0.1015*x.^4);

yc = zeros(size(x)); dyc = yc;
for i=1:length(x)
    if x(i)<p
        yc(i)=m/p^2*(2*p*x(i)-x(i)^2);
        dyc(i)=2*m/p^2*(p-x(i));
    else
        yc(i)=m/(1-p)^2*((1-2*p)+2*p*x(i)-x(i)^2);
        dyc(i)=2*m/(1-p)^2*(p-x(i));
    end
end

thc = atan(dyc);
zu = (x - yt.*sin(thc) )*c_blade;
yu = (yc + yt.*cos(thc))*c_blade;
zl = (x + yt.*sin(thc) )*c_blade;
yl = (yc - yt.*cos(thc))*c_blade;

Zaf = [zu fliplr(zl)];
Yaf = [yu fliplr(yl)];

figure('Name', 'NACA4412') 
hold on
plot(A(:,1), A(:,2), 'Color', 'b', 'DisplayName', 'Point Style')
plot(Zaf, Yaf, 'Color', 'r', 'DisplayName', 'Equation Style')
legend('show')
axis equal
title('NACA 4412 Airfoil');
xlabel('Chord Position');
ylabel('Thickness Distribution');
grid on;
hold off