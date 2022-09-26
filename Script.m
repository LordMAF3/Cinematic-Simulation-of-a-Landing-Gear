load vel_carrello.mat
l1 = 1.2;   % lunghezza dell'asta
l2 = 0.6;   % lunghezza del pistone
l3 = 1.0;   % lunghezza del carrello
d1 = 1.0;   % distanza tra la c.r. di 1 e la c.r. di 3
d2 = 0.15;   % distanza tra la c.r. di 1 e il pistone
v2 = -vel_carrello; % velocità del pistone

T = 3.5;      % durata simulazione
step = 100; % passi della simulazione
vect_scale = 10;

t = linspace(0,T,step)';

% Moto del corpo 2
corsa_pistone = v2.*t;

% Moto del corpo 1
theta = atan(corsa_pistone/d2);

% Moto del corpo 3
BD = [l1*cos(theta)-d1, l1*sin(theta)];
nBD = sqrt(BD(:, 1).^2+BD(:,2).^2);
gamma = asin((l1^2-nBD.^2-d1^2)./(2*d1.*nBD));

%Fine corsa
gamma(gamma<=0)=0;
finecorsa=find(gamma<=0, 1);
theta(finecorsa:step)=theta(finecorsa-1);
corsa_pistone(finecorsa:step)=corsa_pistone(finecorsa-1);
BD(finecorsa:step, 1)=BD(finecorsa-1, 1);
BD(finecorsa:step, 2)=BD(finecorsa-1, 2);

w1 = v2.*cos(theta).^2/d2;
v21 = v2.*sin(theta);
w3 = w1./(1+(d1/l1*(sin(gamma)./sin(theta-gamma))));
v13 = (l1*sin(theta)./sin(gamma)).*(w3-w1);
v13(finecorsa:step) = v13(finecorsa-1);

v_C = [-w1.*d2.*tan(theta), w1.*d2];
v_D = [-w3.*BD(:,2), w3.*BD(:,1)];

%Plot
figure('Name','Cinematica di un carrello');
title('Cinematica di un carrello');
for i=1:step
    axis equal;
    set(gca, 'xlim', [-0.5 2.2], 'ylim',[-1.5 0.7], 'nextplot', 'replacechildren');
    plot([-0.5, 2.2], [0, 0], 'k-', 'LineWidth', 5);
    hold on;
    plot(d1+l3*sin(gamma(i)), -l3*cos(gamma(i)), 'k-o', 'LineWidth', 10, 'MarkerSize', 30); %ruota
    plot([d1 d1+l3*sin(gamma(i))], [0 -l3*cos(gamma(i))], 'r-o', 'LineWidth', 10, 'MarkerSize',5); %carrello
    hold on;
    plot([0 l1*cos(theta(i))],[0 l1*sin(theta(i))],'g-o', 'LineWidth', 5, 'MarkerSize',3); %asta
    hold on;
    plot([d2 d2],[l2+corsa_pistone(i) corsa_pistone(i)],'b', 'LineWidth', 10); %pistone
    plot(d2, corsa_pistone(i), 'b-o','MarkerSize', 5, 'LineWidth', 5)
    hold off;
    hold on;
    quiver(d2, d2*tan(theta(i)), vect_scale*v_C(i,1), vect_scale*v_C(i,2));
    quiver(d2, corsa_pistone(i), 0, vect_scale*v2(i));
    quiver(d2+0.9*vect_scale*v_C(i,1), d2*tan(theta(i))+0.9*vect_scale*v_C(i,2), vect_scale*v21(i).*cos(theta(i)), vect_scale*v21(i).*sin(theta(i)),'k');
    quiver(d1+BD(i,1), BD(i,2), vect_scale/2*v_D(i, 1), vect_scale/2*v_D(i,2));
    quiver(d1+BD(i,1), BD(i,2), vect_scale/2*-w1(i)*l1*sin(theta(i)), vect_scale/2*w1(i)*l1*cos(theta(i)));
    quiver(d1+BD(i,1)+0.9*vect_scale/2*v_D(i, 1), BD(i,2)+0.9*vect_scale/2*v_D(i,2), vect_scale/2*v13(i)*sin(gamma(i)), vect_scale/2*-v13(i)*cos(gamma(i)));
    hold off;
    hold on;
    a=0.11;
    plot([d2-a, d2-a], [0.1, 0.7], 'k', 'LineWidth', 10); %cilindro1
    plot([d2-a, d2+a], [0.655, 0.655], 'k', 'LineWidth', 10); %cilindro2
    plot([d2+a, d2+a], [0.1, 0.7], 'k', 'LineWidth', 10); %cilindro3
    legend('', '', '', '', '', '', 'v_{C-1}', 'v_2', 'v_{21(rel)}', 'v_{D-3}', 'v_{D-1}', 'v_{13(rel)}', '', '', '');
    text(-a, a, 'A');
    text(d1-a, a, 'B');
    text(d1+BD(i,1)-a, BD(i,2)-a, 'D');
    text(d2-a, d2*tan(theta(i))-a, 'C');
    time=sprintf('time=%.3f',t(i));
    text(1.7, -1.4, time);
    hold off;
    try
        F(i)=getframe(gca);
    catch
        break;
    end
end
