function I = linecirc2(ll, cc)
% Find intersection point between all line-circle pairs.
% In the ll matrix each row is a line with format [x1 y1 x2 y2].
% In the cc matrix each row is a circle with format [x y r].

nl = size(ll,1); % Number of lines
nc = size(cc,1); % Number of circles
%nl*nc

% llcc = [repelem(ll, nc, 1) repmat(cc, nl, 1)]; % All combinations
% 
% end1diff = llcc(:,1:2) - llcc(:,5:6);
% end1distance = sqrt(dot(end1diff, end1diff, 2));
% end2diff = llcc(:,3:4) - llcc(:,5:6);
% end2distance = sqrt(dot(end2diff, end2diff, 2));
% 
% llcc(:,8:9) = [(end1distance < llcc(:,7)) (end2distance < llcc(:,7))]

slope = diff( ll(:,[2 4]), 1, 2 ) ./ diff( ll(:,[1 3]), 1, 2 );
intercpt = ll(:,2) - ll(:,1) .* slope;

llcc = [repelem([ll slope intercpt], nc, 1) repmat(cc, nl, 1)]; % All combinations

intersection = nan(nl*nc, 4);
for lc = 1:nl*nc
    [xout, yout] = linecirc( llcc(lc,5), llcc(lc,6), llcc(lc,7), llcc(lc,8), llcc(lc,9) );
    intersection(lc,:) = [xout(1) yout(1) xout(2) yout(2)];
end

% remove nan
I = reshape(intersection(~any(isnan(intersection),2),:)', 2, [])';

    
end