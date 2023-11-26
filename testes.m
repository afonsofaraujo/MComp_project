% Malha de Teste
close all;
clear;
nodesize = 20;

XN = [];

% Block 1
width = 1;
height = 0.9;
xSpacing = 0.1;
ySpacing = 0.1;
[x, y] = meshgrid(0:xSpacing:(width), 0:ySpacing:(height));
XN = [XN; x(:), y(:)];

% Block 2
a = 0.3;  % Semi-major axis
b = 0.6;  % Semi-minor axis
theta = linspace(0, pi, 10);
x = a * cos(theta)+1.3;
y = b * sin(theta);
XN = [XN; x', y']; % Add the new coordinates to XN
numEllipses = 8;
for i = 1:numEllipses
    b = b*0.85 ; % Increase the semi-minor axis
    y = b * sin(theta) + i * 0.1; % Adjust the y-coordinate for each ellipse
        newX = a * cos(theta) + 1.3; % Calculate the x-coordinates for the new ellipse
    XN = [XN; newX', y']; % Add the new coordinates to XN
end
yb2 = 1:xSpacing:1.6;
XN = [XN; yb2(:), 0.9*ones(length(yb2),1)];

% Block 3
width = 1;
height = 0.9;
xSpacing = 0.1;
ySpacing = 0.1;
[x, y] = meshgrid(1.6:xSpacing:(width+1.6), 0:ySpacing:(height));
XN = [XN; x(:), y(:)];

% Block 4
a = 0.3;  % Semi-major axis
b = 0.6;  % Semi-minor axis
theta = linspace(pi/2, pi, 5);
x = a * cos(theta) + 2.9;
y = b * sin(theta);
XN = [XN; x', y']; % Add the new coordinates to XN
numEllipses = 8;
for i = 1:numEllipses
    b = b*0.85 ; % Increase the semi-minor axis
    y = b * sin(theta) + i * 0.1; % Adjust the y-coordinate for each ellipse
    newX = a * cos(theta) + 2.9; % Calculate the x-coordinates for the new ellipse
    XN = [XN; newX', y']; % Add the new coordinates to XN
end
yb4 = 2.6:xSpacing:3;
XN = [XN; yb4(:), 0.9*ones(length(yb4),1)];

% Block 5
width = 1;
height = 0.3;
x0 = 2.9 + xSpacing;
y0 = 0.6;
xSpacing = 0.1;
ySpacing = 0.1;
[x, y] = meshgrid(x0:xSpacing:(width+x0), y0:ySpacing:(height+y0));
XN = [XN; x(:), y(:)];

% Remove duplicates
numPoints = size(XN, 1);
toDelete = false(numPoints, 1);
for i = 1:numPoints
    if toDelete(i)
        continue;
    end
    for j = i+1:numPoints
        tolerance = 1e-6;
        if norm(XN(i, :) - XN(j, :)) < tolerance
            toDelete(j) = true;
        end
    end
end
XN = XN(~toDelete, :);

% Remove points with height greater than 0.9
XN = XN(XN(:, 2) <= 0.9, :);
[Nnds,~] = size(XN);

% % Plot Nodes
% figure;
% scatter(XN(:, 1), XN(:, 2), nodesize, 'filled');
% hold on;
% for L = 1:size(XN, 1)
%     txt = horzcat(' ', num2str(L));
%     text(XN(L, 1) + xSpacing/10, XN(L, 2) - ySpacing/10, txt, 'FontSize', 10);
% end
% hold off;

% Find boundary
k = boundary(XN(:, 1), XN(:, 2),0.73);
numEdges = length(k) - 2;
C = zeros(numEdges, 2);
for H = 1:numEdges
    C(H, :) = [k(H), k(H+1)];
end

DT = delaunayTriangulation(XN);

% Remove wrong triangles
trianglestoremove = [198 199 207 200 232 209 234 415 412 56, 556 567 537 161 553 546 661 663 41 538 557 559 665 674 675]';
trianglesToKeep = setdiff(1:size(DT.ConnectivityList, 1), trianglestoremove);
DT = triangulation(DT.ConnectivityList(trianglesToKeep, :), DT.Points);

% Plot the updated triangulation
figure;
triplot(DT);
axis equal;
for i = 1:size(DT.ConnectivityList, 1)
    centroid = mean(XN(DT.ConnectivityList(i, :), :));
    txt = num2str(i);
    text(centroid(1), centroid(2), txt, 'Color', 'k', 'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
end
Nelt = size(DT,1); % numero de triangulos criados



