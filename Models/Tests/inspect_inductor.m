


mean(cf_0(:, 14:19, :))

format SHORT
v = cf_0(:, 14:15, :);
max(v)
d = sqrt(sum(( v ).^2,2));
format LONG G
max(d)


xy_11_target = [-0.5356 0.4155]; %[1,1]
values = [0.0164 0.0146 2.8994 -0.1834 -0.1493 2.9968];
xy_11_forecast = intercepts + [sum(xy_i_target.*coeffecients(:,[1 3]), 2) sum(xy_i_target.*coeffecients(:,[2 4]), 2)];


xy_i = xy(:,:,i);
xy_past = xy(:,:,i-1);
active_cf = active_cf_i;
cf = cf_i;
a_a = pref.a_a;





unchangedCF = nan(pref.iterations,pref.N);
unactiveCF = nan(pref.iterations,pref.N);
for i=1:pref.iterations
    for n=1:pref.N
        testA = cf(:,14:19,n,1);
        testi = i*pref.psi;
        testB = cf(:,14:19,n,testi);
        %size(testA)
        %size(testB)
        [testIntersect, ~, testib]  = intersect(testA, testB, 'rows');
        textidx = cf(testib,26,n,testi)==0;
        unchangedCF(i, n) = sum(textidx);
        unactiveCF(i, n) = sum( cf(:,26,n,testi)==0 );
    end
end   
unchangedCF
unactiveCF

%diff(unactiveCF)

[testaccuracy, testorder] = sort(cf_i(2:end,24,:), 'descend');
sum(testaccuracy==0)
