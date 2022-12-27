function [alpha, beta, mapping] = phi_construct(k_max, q)

Lambda = 1/sqrt(2);

mapping_pre = containers.Map({0, 1}, {0, 1});
mapping = containers.Map({0, 1}, {0, 1});
for k = 1:k_max
    alpha = zeros(1, 10^(k-1)+1);
    beta = zeros(1, 10^(k-1)+1);
    for i = 1:10^(k-1)+1
        flag = false;
        alpha(i) = (i-1)*10^(-k+1) + 10^(-k) - 2*q*10^(-k);
        beta(i) = (i-1)*10^(-k+1) + 9*10^(-k) - 2*q*10^(-k);
        if alpha(i)>1 || beta(i)<0
            continue
        end
        
        sorted_key = [];
        for gamma = keys(mapping_pre)
            sorted_key = [sorted_key, gamma{1}];
            if alpha(i) < gamma{1} && beta(i) > gamma{1}
                flag = true;
                mapping(alpha(i)) = mapping_pre(gamma{1});
                mapping(beta(i)) = mapping_pre(gamma{1});
            end
        end
        
        %sorted_key = sort(keys(mapping_pre));
        sorted_key = sort(sorted_key);
        size = length(sorted_key);
            
        if ~flag
            for idx_1 = 1:size
                if sorted_key(idx_1) > alpha(i)
                    break
                end
            end
            idx_1 = idx_1 - 1;
            
            for idx_2 = 1:size
                if sorted_key(idx_2) > beta(i)
                    break
                end
            end
            gamma_1 = sorted_key(idx_1);
            gamma_2 = sorted_key(idx_2);
            
            if q == 0
                quotient = floor((alpha(i)-gamma_1)/(10*10^(-k)));
                mapping(alpha(i)) = mapping_pre(gamma_1) + (alpha(i)-gamma_1-quotient*8*10^(-k))*5*(mapping_pre(gamma_2)-mapping_pre(gamma_1))/(gamma_2-gamma_1);
                mapping(beta(i)) = mapping(alpha(i));
            else
                quotient = floor((alpha(i)-gamma_1)/(10*10^(-k)));
                mapping(alpha(i)) = mapping_pre(gamma_1) + (alpha(i)-gamma_1-(9-2*q)*10^(-k)-quotient*8*10^(-k))*5*(mapping_pre(gamma_2)-mapping_pre(gamma_1))/(gamma_2-gamma_1);
                mapping(beta(i)) = mapping(alpha(i));
            end
        end
    end
    
    if q == 0
        if k == 1
            %mapping(alpha(1)) = 13249661/28937850;
            %mapping(beta(1)) = 13373161/28937850;
            mapping(alpha(1)) = 13249661/28937850+0.04;
            mapping(beta(1)) = 13373161/28937850+0.04;
        end
        if k>=2
            d = [];
            %diff_d = [];
            for i = 1:length(alpha)
                for j = 1:length(alpha)
                    if alpha(i) > 1 || alpha(j) > 1
                        continue
                    end
                    d = [d, mapping(alpha(i))+Lambda*mapping(alpha(j))];
                end
            end
            d = sort(d);
            diff_d = zeros(1,length(d)-1);
            for i = 1:length(d)-1
                diff_d(i) = d(i+1)-d(i);
                %diff_d = [diff_d, d(i+1)-d(i)];
            end
            
            epsilon = min(diff_d)/8;
            %epsilon = 0
            for i = 1:length(alpha)
                if alpha(i) > 1 || beta(i) < 0
                    continue
                end
                mapping(alpha(i)) = mapping((alpha(i))) - epsilon;
                mapping(beta(i)) = mapping(beta(i)) + epsilon;
            end
        end
    else
        if k == 1
            %mapping(alpha(1)) = -(13373161/28937850-0.46);
            %mapping(beta(1)) = 13373161/28937850-0.46;
            %mapping(alpha(2)) = 23/25 - (13373161/28937850-0.46);
            %mapping(beta(2)) = 23/25 + (13373161/28937850-0.46);
            mapping(alpha(1)) = -(13373161/28937850-0.46);
            mapping(beta(1)) = 13373161/28937850-0.46;
            mapping(alpha(2)) = 1 - (13373161/28937850-0.46);
            mapping(beta(2)) = 1 + (13373161/28937850-0.46);
        end
        if k>=2
            d = [];
            %diff_d = [];
            for i = 1:length(alpha)
                for j = 1:length(alpha)
                    if alpha(i) > 1 || alpha(j) > 1
                        continue
                    end
                    d = [d, mapping(alpha(i))+Lambda*mapping(alpha(j))];
                end
            end
            d = sort(d);
            diff_d = zeros(1, length(d)-1);
            for i = 1:length(d)-1
                diff_d(i) = d(i+1)-d(i);
                %diff_d = [diff_d, d(i+1)-d(i)];
            end

            
            epsilon = min(diff_d)/8;
            %epsilon
            %epsilon = 0
            for i = 1:length(alpha)
                if alpha(i) > 1 || beta(i) < 0
                    continue
                end
                mapping(alpha(i)) = mapping(alpha(i)) - epsilon;
                mapping(beta(i)) = mapping(beta(i)) + epsilon;
            end
        end
    end
    %length(mapping_pre)
    %length(mapping)
    
    mapping_pre = containers.Map(mapping.keys, mapping.values);
    
end 

end

