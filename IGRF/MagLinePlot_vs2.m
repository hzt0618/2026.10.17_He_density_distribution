%% ÅäÌ×MagLine_vs2
figure(1);

load line_glat_30.mat

for i=1:size(line,2)
plot3(line{i}(:,1),line{i}(:,2),line{i}(:,3));
hold on;
end

[X,Y,Z] = sphere(100);
r = 6371.2;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

s=surf(X2,Y2,Z2);
s.EdgeColor = 'none';
axis equal
hold on