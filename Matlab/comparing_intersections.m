p1 = randi(15,3,1);
p2 = randi(15,3,1);
p3 = randi(15,3,1);
points = [p1,p2,p3];

y = 1;
z = 1;

normal = cross(points(:,2)-points(:,1),points(:,3)-points(:,1));
d = dot(normal,points(:,3));
if normal(1) ~= 0
    x = (d-(normal(2)*y)-(normal(3)*z))/normal(1);
else
    fprintf("no intersection")
end

t = (z-p3(3)-(((y-p2(2))*(p3(3)-p2(3)))/(p3(2)-p2(2))))/((p1(3)-p2(3))-(((p1(2)-p2(2))*(p3(3)-p2(3)))/(p3(2)-p2(2))));
s = ((y-p2(2))/(p3(2)-p2(2)))-t*((p1(2)-p2(2))/(p3(2)-p2(2)));
x_prime = p2(1) + s*(p3(1)-p2(1)) + t*(p1(1)-p2(1));

fprintf("normal vector method")
x
p = [x;y;z];
fprintf("linear combo method")
x_prime
p_prime = [x_prime,y,z];

figure;
hold on;
fill3([p1(1), p2(1), p3(1)], [p1(2), p2(2), p3(2)], [p1(3), p2(3), p3(3)], 'b', 'FaceAlpha', 0.3);
scatter3(p(1), p(2), p(3), 100, 'r', 'filled');
scatter3(p_prime(1), p_prime(2), p_prime(3), 100, 'g', 'filled');
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Triangle and Points');
axis equal;
grid on;
legend('Plane Simplex', 'Point p', 'Point p_prime');
hold off;
