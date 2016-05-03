

points = [-0.1325 -2.2267; -0.1525 -2.2267; -0.5319  1.0698; -1.3628 -0.1296;  1.7438  1.3784;  1.5770  0.9458;  0.5147 -2.6114;  0.8169 -2.2797; -1.0244  2.7143; -0.4422  2.8257; -1.7421 -2.4453; -2.4492 -0.4012]
linesegments = [-1.1258 -4.2270 -0.7196 -3.9662; -0.7196 -3.9662  0.4347 -0.4873; -2.3293  1.4275 -3.3717  2.2654; -2.3293  1.4275  0.4347 -0.4873;  1.3579  3.1700  3.3566  0.5079;  3.3566  0.5079  0.4347 -0.4873] % Each row is line with format [x1 y1 x2 y2];

ls = reshape([linesegments nan(size(linesegments,1),2)]', 2, [])';
ch = convhull( points(:,1), points(:,2) );



figure(15);
plot( points(ch,1) , points(ch,2), 'Color', [1 0.64 0]);
hold on;
scatter( points(:,1), points(:,2), [], 'k', 'filled');
text( points(:,1), points(:,2)+0.25, cellstr(num2str([1:size(points,1)]')), 'Color','black');
scatter( ls(:,1), ls(:,2) , 'b');
plot( ls(:,1), ls(:,2) , 'b');
xlim([-5 5]); ylim([-5 5]);
hold off;

ch1 = [9 10 5 6 3 9];
ch2 = [12 4 2 11 12];
ch3 = [1 8 7 1];


figure(16);
plot( points(ch1,1) , points(ch1,2), 'Color', [1 0.64 0]);
hold on;
plot( points(ch2,1) , points(ch2,2), 'Color', [1 0.64 0]);
plot( points(ch3,1) , points(ch3,2), 'Color', [1 0.64 0]);
scatter( points(:,1), points(:,2), [], 'k', 'filled');
text( points(:,1), points(:,2)+0.25, cellstr(num2str([1:size(points,1)]')), 'Color','black');
scatter( ls(:,1), ls(:,2) , 'b');
plot( ls(:,1), ls(:,2) , 'b');
xlim([-5 5]); ylim([-5 5]);
hold off;