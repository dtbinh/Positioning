% Get all the convex sets with at least two triangles.
[L, I] = sort(cellfun(@length,convexTriangleSets),'descend');
L_idx = I(L>1);

combi = cell(0);
% Loop through these convex sets.
for ll = 1:length(L_idx)
    % Extract array with triangles
    A = convexTriangleSets{L_idx(ll)};
    % Calculate all using these triangles. (don't include full array itself).
    for i = 1:size(A,2)-1
        combi = [combi; num2cell(nchoosek(A,i),2)];
    end
end
% Remove duplicates
[~,c] = unique(cellfun(@(x) sprintf('%d,',x),combi,'UniformOutput',false))
combi = combi(c);

% Remove the convex sets that are not the largest.
rm = ismember(cellfun(@(x) sprintf('%d,',x),convexTriangleSets','UniformOutput',false), cellfun(@(x) sprintf('%d,',x),combi,'UniformOutput',false));
convexTriangleSets(rm) = [];


