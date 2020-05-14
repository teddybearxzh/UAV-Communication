%%
close all
clear all
clc

K = 3;     %cell number
M = 4;     %users number
r = 500;   %cell radius

%parameters
n_dis = 5;
gamma = 4;
sigma_square = power(10, -107/10);
P_num = 10;

%cp1 = [0, r];
cp1 = [0, r];
cp2 = [-sqrt(3)*r/2, - r/2];
cp3 = [sqrt(3)*r/2, - r/2];

cp_total = {cp1, cp2, cp3};
DAE_total = cell(1, K);
MT_total = cell(1, K);

for i = 1 : K
    [DAE, MT, d] = AE_MT_generating(cp_total, i, M, r, n_dis);
    DAE_total{i} = DAE;
    MT_total{i} = MT;
end

%%
figure;
hold on;
for i = 1 : K
	%scatter(cp_total{i}(1), cp_total{i}(2), 'r', 's');
	f = ['(x-', num2str(cp_total{i}(1)), ')^2+(y-', num2str(cp_total{i}(2)),...
        ')^2-500^2'];
    ezplot (f, [-1200,1200]);
    for j = 1 : M
    	scatter(DAE_total{i}{j}(1), DAE_total{i}{j}(2), 'b');
    end
    scatter(MT_total{i}(1), MT_total{i}(2), 'm');
end
axis equal;
title('Topology');

total_cell = cell(1, K);
for i = 1 : K
    new_cell = cell(1, 3);
    new_cell{1} = cp_total{i};
    new_cell{2} = DAE_total{i};
    new_cell{3} = MT_total{i};
    total_cell{i} = new_cell;
end
for i = 1 : K
    for j = 1 : M
    	scatter(total_cell{i}{2}{j}(1), total_cell{i}{2}{j}(2), 'b');
        if j == M
            plot([total_cell{i}{2}{j}(1),total_cell{i}{2}{1}(1)],[total_cell{i}{2}{j}(2),total_cell{i}{2}{1}(2)], 'b--');
        else
            plot([total_cell{i}{2}{j}(1),total_cell{i}{2}{j+1}(1)],[total_cell{i}{2}{j}(2),total_cell{i}{2}{j+1}(2)], 'b--');
    
        end
    end
end
clear cp1 cp2 cp3
clear MT DAE new_cell
clear f i j

%%
H = 150;
HH = zeros(K,M);
for i = 1 : K
    for j = 1 : M
        HH(i,j) = sqrt(g(total_cell, i, i, j, H, M));
    end
end
%H and L
% H = cell(K, K);
% L = cell(K, K);
% for k = 1 : K
%     for i = 1 : K
%         %H_ki_w
%         H_ki_w = zeros(M, M);
%         for m = 1 : M
%             for n = 1 : M
%                 H_ki_w(m, n) = complex(normrnd(0,1), normrnd(0,1)) / sqrt(2);
%             end
%         end
%         %L_ki
%         L_ki = zeros(M, M);
%         for n = 1 : M
%             L_ki(n, n) = power(norm(total_cell{i}{2}{n} - total_cell{k}{3}), -gamma) * power(10, normrnd(0, 8) / 10);
%         end
%         L{k, i} = L_ki;
%         H{k, i} = H_ki_w * L_ki;
%     end
% end

%%
E_total = [];
data = zeros(3,P_num);
for pp = 1 : P_num
    E = 1 + 99 / (P_num - 1) * (pp - 1);
    E_total = [E_total, E];
    %E = power(10, E / 10);
    
    
    %W1
    p1 = E / K * ones(1, M);
    p_average = E / K * ones(K, M);
    if E/K >= 1000
        p_average = 1000 * ones(K, M);
    end
    W1 = 1;
    e = 0.000001;
    while(1)
        W_new1 = 0;
        for n = 1 : M
            W_new1 = W_new1 + (HH(1,n)^2*p1(n))/(sigma_square + HH(1,n)^2*p1(n)/W1*M);
        end
        W_old1 = W1;
        W1 = W_new1 + 1;
        if((W1-W_old1)^2) < e
            break;
        end
    end
    
    %W2
    p2 = E / K * ones(1, M);
    W2 = 1;
    e = 0.000001;
    d = [0.95,0.78,0.81,0.86,0.90,0.93,0.95,0.97,0.99,0.998];
    while(1)
        W_new2 = 0;
        for n = 1 : M
            W_new2 = W_new2 + (HH(2,n)^2*p2(n))/(sigma_square + HH(2,n)^2*p2(n)/W2*M);
        end
        W_old2 = W2;
        W2 = W_new2 + 1;
        if((W2-W_old2)^2) < e
            break;
        end
    end
    
    %W3
    p3 = E / K * ones(1, M);
    W3 = 1;
    e = 0.000001;
    while(1)
        
        W_new3 = 0;
        for n = 1 : M
            W_new3 = W_new3 + (HH(3,n)^2*p3(n))/(sigma_square + HH(3,n)^2*p3(n)/W3*M);
        end
        W_old3 = W3;
        W3 = W_new3 + 1;
        if((W3-W_old3)^2) < e
            break;
        end
    end
    
    s = 1;
    while(1)
        s = s + 1;
        %update p
        cvx_begin
            variables p(K,M)
            C1 = 0;
            C1_average = 0;
            for i = 1 : M
                C1 = C1 + log(1+HH(1,i)^2*p(1,i)/W1*M/sigma_square);
                C1_average = C1_average + log(1+HH(1,i)^2*p_average(1,i)/W1*M/sigma_square);
            end
            C1 = C1 + M*log(W1) - M*log(exp(1)*(1-1/W1));
            C1_average = C1_average + M*log(W1) - M*log(exp(1)*(1-1/W1));
            
            C2 = 0;
            C2_average = 0;
            for i = 1 : M
                C2 = C2 + log(1+HH(2,i)^2*p(2,i)/W2*M/sigma_square);
                C2_average = C2_average + log(1+HH(2,i)^2*p_average(2,i)/W2*M/sigma_square);
            end
            C2 = C2 + M*log(W2) - M*log(exp(1)*(1-1/W2));
            C2_average = C2_average + M*log(W2) - M*log(exp(1)*(1-1/W2));
            
            C3 = 0;
            C3_average = 0;
            for i = 1 : M
                C3 = C3 + log(1+HH(3,i)^2*p(3,i)/W3*M/sigma_square);
                C3_average = C3_average + log(1+HH(3,i)^2*p_average(3,i)/W3*M/sigma_square);
            end
            C3 = C3 + M*log(W3) - M*log(exp(1)*(1-1/W3));
            C3_average = C3_average + M*log(W3) - M*log(exp(1)*(1-1/W3));
            
            C=C1+C2+C3;
            C_average = C1_average+C2_average+C3_average;
            
            maximize(C)
            subject to
                sum(p) <= [E,E,E,E]
                for i = 1 : K
                    for j = 1 : M
                        p(i,j) >= 0
                        p(i,j) <= 1000
                    end
                end
        cvx_end
        
        %updateW1W2W3
        %W1
        W1 = 1;
        e = 0.000001;
        %d = [11,11.5,12.7,14.5,15.2,15.5]/11;
        while(1)
            W_new1 = 0;
            for n = 1 : M
                W_new1 = W_new1 + (HH(1,n)^2*p(1,n))/(sigma_square + HH(1,n)^2*p(1,n)/W1*M);
            end
            W_old1 = W1;
            W1 = W_new1 + 1;
            if((W1-W_old1)^2) < e
                break;
            end
        end

        %W2
        W2 = 1;
        e = 0.000001;
        while(1)
            W_new2 = 0;
            for n = 1 : M
                W_new2 = W_new2 + (HH(2,n)^2*p(2,n))/(sigma_square + HH(2,n)^2*p(2,n)/W2*M);
            end
            W_old2 = W2;
            W2 = W_new2 + 1;
            if((W2-W_old2)^2) < e
                break;
            end
        end

        %W3
        W3 = 1;
        e = 0.000001;
        while(1)
            W_new3 = 0;
            for n = 1 : M
                W_new3 = W_new3 + (HH(3,n)^2*p(3,n))/(sigma_square + HH(3,n)^2*p(3,n)/W3*M);
            end
            W_old3 = W3;
            W3 = W_new3 + 1;
            if((W3-W_old3)^2) < e
                break;
            end
        end
        
        if s == 3
            display('*************');
            break;
        end
    end
    Pt_1 = diag(p_average(1,:));
    Pt_2 = diag(p_average(2,:));
    Pt_3 = diag(p_average(3,:));
    data(2,pp) = (log2(det(eye(M)+HH(1,:)*Pt_1*HH(1,:)'/sigma_square)) + log2(det(eye(M)+HH(2,:)*Pt_2*HH(2,:)'/sigma_square)) + log2(det(eye(M)+HH(3,:)*Pt_3*HH(3,:)'/sigma_square)))/2;
    data(1,pp) = C;
    data(3,pp) = C_average;
end
%%
figure;
%Es = power(10, E_total ./ 10)/1000;
%semilogx(E_total,abs(data(1,:)/2));
hold on;
%semilogx(E_total,data(2,:)/2.*d);
ylabel('Total Ergodic capacity(bits/s/Hz)');
xlabel('Transmit Energy constraint for a single UAV (J)');
plot(E_total,abs(data(1,:)),'bo-');
plot(E_total,data(3,:).*d,'rx-');
%plot(E_total,data(2,:)/max(data(2,:))*max(data(1,:)),'gx-');
grid on;
legend('Proposed Scheme','Average Scheme','Location','Best');


