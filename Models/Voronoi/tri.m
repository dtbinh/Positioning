xy_boundary = [3.3735 0.7889; -0.1072 -3.4814; -3.9732 4.1955; -5 5; 5 5; 5 -5; -5 -5];
%xy_boundary = [xy(:,:,end); -5 5; 5 5; 5 -5; -5 -5]

[convexTriangleSets,DT] = triangles(xy_boundary);