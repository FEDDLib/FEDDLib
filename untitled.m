E0 = 1000;
E1 = 900;
c1 = 0.1;
c = linspace(0,1,1000);

f_vec1 = E0 -(E0-E1)*(c./(c+c1));

c1 = 0.01;
f_vec2 = E0 -(E0-E1)*(c./(c+c1));

c1 = 0.001;
f_vec3 = E0 -(E0-E1)*(c./(c+c1));

E0(1:1000) = E0;
E1(1:1000) = E1;
semilogy(c,f_vec1,'g ',...
'LineWidth',1.5,...
'MarkerSize',9, ...
'MarkerEdgeColor','g');
hold on
semilogy(c,f_vec2,'b ',...
'LineWidth',1.5,...
'MarkerSize',9, ...
'MarkerEdgeColor','b');
hold on
semilogy(c,f_vec3,'r ',...
'LineWidth',1.5,...
'MarkerSize',9, ...
'MarkerEdgeColor','r');
hold on
plot(c,E0, '-- k','LineWidth',1.0);
hold on
plot(c,E1,' : k ','LineWidth',1.5);
xlabel('Concentration c')
xlim([-0.01 1])
ylim([895 1005])
ylabel('Youngs Modulus')
legend('c_1=0.1','c_1=0.01','c_1=0.001', 'E_0', 'E_1');
grid on