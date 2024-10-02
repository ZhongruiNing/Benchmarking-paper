function [bestx, calibration_result] = sceua(fn, x0, bl, bu, maxn, ngs, iseed, x_obs, y_obs, fn_hm)
    % Initialize SCE parameters:
    nopt = length(x0); % 参数个数
    npg  = 2 * nopt + 1; % 参数2倍+1，复合体（亚种群数量）
    nps  = nopt + 1; % 单形中的成员数
    nspl = npg;  % 洗牌前每个复合体的进化步数
    npt  = npg * ngs; % 种群数量

    % make sure x is a row matrix
    x0 = as_rowvec(x0);
    bu = as_rowvec(bu);
    bl = as_rowvec(bl);

    bound = bu - bl; % 参数取值范围

    % Create an initial population to fill array x(npt,nopt):
    rand('seed', iseed);

    x = zeros(npt,nopt); % 生成一个长度为npt，宽度为参数个数的全0数组，用来存放随机生成的参数

    for i = 1 : npt
        x(i, :) = bl + rand(1, nopt) .* bound; % 随机生成参数
    end

    x(1,:) = x0;


    for i = 1 : npt
        xf(i) = fn(x_obs, x(i,:), y_obs, fn_hm); % nopt，根据随机生成的参数计算损失函数
    end
    f0=xf(1); % 取第一个损失函数值作为初始值

    % Sort the population in order of increasing function values;
    [xf, idx] = sort(xf); % 对损失函数值从小到大排序，并获取顺序序列
    x = x(idx, :); % 根据损失函数值的排序结果，对函数数组进行排序。

    % Record the best and worst points;
    bestx  = x(1,:); 
    bestf  = xf(1); %记录最好的参数数值与损失函数值
    BESTF  = bestf;
    BESTX  = bestx;

    % Compute the standard deviation for each parameter
    % xnstd = std(x); % 计算每个参数的标准差（npt个数值计算标准差）

    % Computes the normalized geometric range of the parameters
    % gnrng = exp(mean(log((max(x) - min(x)) ./ bound))); % 计算参数的归一化几何范围


    % Begin evolution loops: 开始进化循环
    % criter = [];

    % 循环的条件，icall小于最大试验次数maxn，参数范围gnrng大于预先设定的peps，criter_change大于pcento
    cr_n = 1;
    for tt = 1 : maxn
        % Loop on complexes (sub-populations); 复合体上的循环，亚种群
        % 将种群分为ngs个复形，每个复形中有npg个单形。每次循环一个复形
        % 每个复形进化nspl次，每次改变最差的一个单形，将进化后的单形重新混合进复形
        for igs = 1 : ngs
            % Partition the population into complexes (sub-populations); 将种群划分为复合体(亚种群);
            k1 = 1 : npg; 
            k2 = (k1 - 1) * ngs + igs; % 将所有种群分为ngs个复合体，每个复合体里面又npg个亚种
            cx(k1,:) = x(k2,:); % 复合体内所有亚种的参数
            cf(k1)   = xf(k2); % 复合体内所有亚种的损失函数值
            % Evolve sub-population igs for nspl steps: 进化亚种群igs的nspl步骤:
            for loop = 1 : nspl % 循环nspl次（进化nspl次），每次进化都会改变cx中最差的一个点
                % Select simplex by sampling the complex according to a linear probability distribution
                % 根据线性概率分布对复形进行抽样，选择单纯形
                % 从1到npg+1取出nps个不重复整数，存储在lcs中
                lcs(1) = 1;
                for k3=2:nps
                    for iter=1:1000
                        % 取一个从1到npg的随机整数
                        lpos = 1 + floor(npg+0.5-sqrt((npg+0.5)^2 - npg*(npg+1)*rand));
                        % 如果这个数已经存在于lcs中，则重新取，否则，这个数就是需要取出的数
                        idx=find(lcs(1:k3-1)==lpos, 1); 
                        if isempty(idx); break; end
                    end
                    lcs(k3) = lpos;
                end
                % 对lcs排序
                lcs=sort(lcs);

                % Construct the simplex: 构建单纯形
                s = zeros(nps,nopt);
                % 从复合体内取出单纯性位置的亚种，并获取他们的损失函数值
                s=cx(lcs,:); sf = cf(lcs);
                % cce算法，算法的作用是根据取出的亚种，计算除了最坏的点之外的点的质心，之后根据质心
                % 计算一个收缩点和一个反射点，判断哪个点更好，就用哪个点了，如果都更差，则生成一个随机的点
                [snew, fnew] = cceua(fn, s, sf, bl, bu, x_obs, y_obs, fn_hm);

                % Replace the worst point in Simplex with the new point:
                % 用cce算法中得到的新点替换单形中最差的点
                s(nps,:) = snew; sf(nps) = fnew;

                % Replace the simplex into the complex;
                % 将更新后的单形还原到复形中，实际上只改变了最后一个点
                cx(lcs,:) = s;
                cf(lcs) = sf;

                % Sort the complex;
                % 重新排列复形
                [cf,idx] = sort(cf); cx=cx(idx,:);
                % End of Inner Loop for Competitive Evolution of Simplexes
                % 单纯形竞争进化的内循环结束
                calibration_result(cr_n) = xf(1);
                cr_n = cr_n + 1;
            end
            % Replace the complex back into the population;
            % 更新后的复形重新纳入种群
            x(k2,:) = cx(k1,:);
            xf(k2) = cf(k1);
            % End of Loop on Complex Evolution; 结束复形进化循环
            
            
        end
        % Shuffled the complexes; 复合体洗牌
        % 对进化后的种群进行排序，并将排序后的种群及其对应的损失函数值赋值给PX和PF
        [xf,idx] = sort(xf); x=x(idx,:);


        % Record the best and worst points; 记录下最好的点和最坏的点
        bestx=x(1,:); bestf=xf(1);

        % 和初始阶段最好的点和最坏的点进行混合
        BESTX=[BESTX;bestx]; BESTF=[BESTF;bestf];

        % Compute the standard deviation for each parameter
        % 计算进化后种群每个参数的标准差
        % xnstd=std(x);

        % Computes the normalized geometric range of the parameters
        % 计算进化后种群参数的质心
        gnrng=exp(mean(log((max(x)-min(x))./bound)));

        % criter=[criter;bestf];
        % End of the Outer Loops
        if gnrng < 0.0001
            return
        end

        tt = length(calibration_result);
        if tt > 1
            if (calibration_result(tt) - calibration_result(tt-1)) < 0.01
                return
            end
        end
    end

end

function bool = is_rowvec(x)
    [nrow, ncol] = size(x);
    bool = nrow == 1 && ncol > 1;
end

function bool = is_colvec(x)
    [nrow, ncol] = size(x);
    bool = ncol == 1 && nrow > 1;
end

function x = as_rowvec(x)
    if is_rowvec(x)
       % nothing 
    elseif is_colvec(x)
        x = x';  
    end
end

function x = as_colvec(x)
    if is_colvec(x)
       % nothing 
    elseif is_rowvec(x)
        x = x';  
    end
end