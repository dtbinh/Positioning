
[L,I] = sort(cellfun(@length,convexTriangleSets),'descend');
[~, L_idx] = min(L)

convexTriangleSets{I(1:L_idx-1)}



[~, adj_break] = max(cumsum(adj,1), [], 1);
adj_break = [0 adj_break];

for tri = 1:size(adj,2)
    adj_tri = adj(adj_break(tri)+1:adj_break(tri+1),:)
    pdist(adj_tri)
    
end




test = repmat(adj, size(adj,1), 1)-repelem(adj, size(adj,1), 1)

sum_adj = sum(adj,2,'omitnan')


% Subtracting the conditions from the state J gives 0 if a particular
% condition is statisfied, and 1 or -1 if it is not satisfied. Summing over
% all the particular conditions reveals if any of the conditions are 
% unfulfilled. By ignoring NaN when summing makes NaN a wildcard character, 
% that will match both 0s and 1s in state J.
count_unfulfilled = sum( abs( repmat(J, M, 1)-cf(:,1:13) ), 2, 'omitnan' );
% Index of conditions that satisfy the current state J.
c_idx = find( count_unfulfilled == 0);







%%%%%% 




teststr = '1,3,6,';
testarray = [1 3 6];

A = testarray;

for i = 2:size(A,2)-1
    combi = [combi; num2cell(nchoosek(A,i),2)];
end







teststring = cellfun(@(x) sprintf('%d,',x), convexTriangleSets, 'UniformOutput',false)

regexp(teststring, '3|6')

