function D = calcProfile(A1,A2,varargin)
% This function computes connectivity profile distance between two matching
% matrices A1 and A2.

    connLimit = [];

    % Parse optional arguments
    while ~isempty(varargin)
        if numel(varargin) == 1
            error('lscatter:missing_option', ...
                'Optional arguments must come in pairs.');
        end

        switch lower(varargin{1})
            case 'connlimit'
                assert(isnumeric(varargin{2}))
                connLimit = varargin{2};                     
        end

        varargin(1:2) = [];

    end


    if size(A1) ~= size(A2)
        error('Dimensions of A1 and A2 must match.')
    end

    N = size(A1, 1);
    D = nan(N);

    for i = 1:N
        
        v1 = A1(i,:);
        
        if ~isempty(connLimit)
            
            [~, I] = sort(v1, 'descend');
            v1(I(connLimit+1:end)) = 0; % only keep the top connLimit strongest connections
            
        end  
        
        for j = 1:N
            
            tmp1 = v1;
            tmp2 = A2(j,:);
            
            if ~isempty(connLimit)

                [~, I] = sort(tmp2, 'descend');
                tmp2(I(connLimit+1:end)) = 0; % only keep the top connLimit strongest connections

            end            

            tmp1([i,j]) = []; % remove self-connections and connections with j
            tmp2([i,j]) = []; % remove self-connections and connections with i

            % Calculate connectivity profile distance
            tmp3=abs(tmp2-tmp1);
            tmp3(tmp3==0) = NaN;
            D(i,j) = nanmean(tmp3);
            
        end

    end
end