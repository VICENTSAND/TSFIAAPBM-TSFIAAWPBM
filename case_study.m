% 数字营销案例：动态支持 IAA--BM 算子。
% 直接点击运行；输入数据、计算过程和结果输出均在本文件中。

%% 1. 输入评价矩阵
caseData = zeros(5,8,3);
% 三个矩阵依次为隶属度 S、非隶属度 D、中立度 H；行是方案，列是准则。
caseData(:,:,1)=[
    0.95 0.85 0.65 0.75 0.95 0.85 0.95 0.65;
    0.75 0.55 0.55 0.75 0.65 0.85 0.65 0.35;
    0.95 0.85 0.65 0.55 0.65 0.95 0.85 0.55;
    0.85 0.65 0.65 0.75 0.65 0.85 0.75 0.55;
    0.75 0.55 0.55 0.65 0.45 0.55 0.65 0.45
];
caseData(:,:,2)=[
    0.05 0.15 0.25 0.20 0.05 0.15 0.05 0.25;
    0.20 0.35 0.35 0.20 0.25 0.15 0.25 0.55;
    0.05 0.15 0.25 0.35 0.25 0.05 0.15 0.35;
    0.15 0.25 0.25 0.20 0.25 0.15 0.20 0.35;
    0.20 0.35 0.35 0.25 0.45 0.35 0.25 0.45
];
caseData(:,:,3)=[
    0.05 0.10 0.20 0.15 0.05 0.10 0.05 0.20;
    0.15 0.30 0.30 0.15 0.20 0.10 0.20 0.50;
    0.05 0.10 0.20 0.30 0.20 0.05 0.10 0.30;
    0.10 0.20 0.20 0.15 0.20 0.10 0.15 0.30;
    0.15 0.30 0.30 0.20 0.40 0.30 0.20 0.40
];

%% 2. 设置参数和两组准则权重
gamma = 3;
alpha = 1;
P = 1;
Q = 1;
tau = 0.5;
kappa = 0.5;
weightSets = [ones(1,8)/8;
              0.20 0.15 0.10 0.10 0.15 0.10 0.10 0.10];
settingNames = {'equal','unequal'};
benefitCriteria = true(1,8);        % 本案例的八个准则均为效益型

tolerance = 5e-13;                  % 仅用于浮点舍入误差检查
assert(all(isfinite([gamma,alpha,P,Q,tau,kappa])) && ...
    gamma >= 1 && alpha >= 1 && P > 0 && Q > 0 && isfinite(P+Q) && ...
    tau >= 0 && tau < 1 && kappa >= 0 && kappa <= 0.5, ...
    'TSF:Parameter','Invalid parameters.');
assert(all(isfinite(caseData),'all') && all(caseData >= 0 & caseData <= 1,'all'), ...
    'TSF:Input','Component degrees must lie in [0,1].');
[nAlternatives,nCriteria,~] = size(caseData);
assert(numel(benefitCriteria) == nCriteria && size(weightSets,2) == nCriteria);

%% 3. 统一准则方向，转换为幂坐标
% 成本型准则交换 S 与 D，中立度 H 不变。
normalizedData = caseData;
normalizedData(:,~benefitCriteria,[1,2]) = caseData(:,~benefitCriteria,[2,1]);
mass = normalizedData.^gamma;
assert(~any(normalizedData > 0 & mass == 0,'all'), ...
    'TSF:Precision','Component-power conversion underflowed.');
assert(all(sum(mass,3) <= 1+tolerance,'all'), ...
    'TSF:Input','Inadmissible T-spherical fuzzy input.');

alternatives = compose('M%d',(1:nAlternatives)');
nSettings = size(weightSets,1);
allAggregates = zeros(nAlternatives,3,nSettings);
allIndices = zeros(nAlternatives,3,nSettings);
caseTables = cell(1,nSettings);
caseRankings = cell(1,nSettings);

%% 4. 分别计算等权和不等权结果
for setting = 1:nSettings
    w = weightSets(setting,:);
    assert(all(isfinite(w)) && all(w > 0) && abs(sum(w)-1) < 2e-14, ...
        'TSF:Weights','Criterion weights must be positive and normalized.');

    % 完整有序对包含 i=j，配对系数为 w(i)*w(j)。
    pairWeights = w.' * w;
    assert(all(pairWeights > 0,'all') && abs(sum(pairWeights,'all')-1) < 2e-14, ...
        'TSF:Precision','Pair coefficients are not representable.');
    scale = (P+Q)^(1/alpha);
    assert(scale > 0 && isfinite(scale), ...
        'TSF:Precision','Pair scale is not representable.');

    aggregate = zeros(nAlternatives,3);
    indices = zeros(nAlternatives,3);

    for a = 1:nAlternatives
        % 4.1 取出一个方案的 s=S^gamma、d=D^gamma、h=H^gamma。
        s = mass(a,:,1).';
        d = mass(a,:,2).';
        h = mass(a,:,3).';
        cumulativeColumns = [s, 1-d, s+h];

        % 4.2 在三个累积坐标上，先求动态支持，再同步交换。
        exchanged = zeros(nCriteria,3);
        for channel = 1:3
            z = cumulativeColumns(:,channel);
            difference = z.' - z;                  % 第 (i,j) 项为 z(j)-z(i)
            support = 1-kappa*abs(difference);
            exchanged(:,channel) = z + tau*((difference.*support)*w.');
        end
        sEffective = exchanged(:,1);
        dEffective = 1-exchanged(:,2);
        hEffective = exchanged(:,3)-exchanged(:,1);
        effective = [sEffective,dEffective,hEffective];
        assert(all(effective >= -tolerance,'all') && ...
            all(sum(effective,2) <= 1+tolerance), ...
            'TSF:Closure','Invalid effective triple.');
        effective = min(max(effective,0),1);

        % 4.3 对每一对准则计算 IAA 生成元。
        % R(i,j)=(P*(-log(z(i)))^alpha+Q*(-log(z(j)))^alpha)^(1/alpha)。
        % 在对数域计算同一公式，避免直接求幂时溢出。
        cumulativeColumns = [effective(:,1), 1-effective(:,2), ...
                             min(effective(:,1)+effective(:,3),1)];
        pairGenerators = zeros(nCriteria,nCriteria,3);
        for channel = 1:3
            z = cumulativeColumns(:,channel);
            g = -log(z);
            left = log(g)+log(P)/alpha;
            right = log(g.')+log(Q)/alpha;
            peak = max(left,right);
            logGenerator = peak + log(exp(alpha*(left-peak)) + ...
                                     exp(alpha*(right-peak)))/alpha;
            logGenerator(isinf(peak)) = peak(isinf(peak));
            generator = exp(logGenerator);
            positivePair = (z > 0) & (z.' > 0);
            assert(~any(isnan(generator) | (isinf(generator) & positivePair) | ...
                (isfinite(generator) & generator > 700),'all'), ...
                'TSF:Precision','Pair generator exceeds the supported full-triple precision range.');
            pairGenerators(:,:,channel) = generator;
        end

        % 4.4 还原配对结果，保留中立度的差值项。
        membershipGenerator = pairGenerators(:,:,1);
        nonmembershipGenerator = pairGenerators(:,:,2);
        cumulativeGenerator = pairGenerators(:,:,3);
        gap = membershipGenerator-cumulativeGenerator;
        gap(isinf(membershipGenerator) & isinf(cumulativeGenerator)) = 0;
        assert(~any(gap < -tolerance.*max(1,abs(membershipGenerator)),'all'), ...
            'TSF:Precision','Negative pair-neutrality difference.');
        % expm1(x)=exp(x)-1；log1p(x)=log(1+x)，小量计算更准确。
        pairNeutrality = exp(-cumulativeGenerator).*(-expm1(-max(gap,0)));
        pairNonmembership = -expm1(-nonmembershipGenerator);
        pairCumulative = pairNonmembership+pairNeutrality;
        assert(all(pairCumulative <= 1+tolerance,'all'), ...
            'TSF:Closure','Invalid pair cumulative component.');
        generators = cat(3, -log(-expm1(-membershipGenerator)), ...
            -log(pairNonmembership), -log(min(pairCumulative,1)));

        % 4.5 按配对权重汇总：N=(sum(w(i)*w(j)*g(i,j)^alpha))^(1/alpha)。
        % 除以最大值后再求幂，公式不变，也避免中间量溢出。
        norms = zeros(1,3);
        for channel = 1:3
            g = generators(:,:,channel);
            peak = max(g,[],'all');
            if isinf(peak)
                norms(channel) = Inf;
            elseif peak == 0
                norms(channel) = 0;
            else
                weightedSum = sum(pairWeights.*(g/peak).^alpha,'all');
                norms(channel) = peak*weightedSum^(1/alpha);
            end
        end

        % 4.6 完成 BM 外层变换，得到聚合后的三个幂坐标。
        gap = norms(2)-norms(3);
        if isinf(norms(2)) && isinf(norms(3))
            gap = 0;
        end
        assert(gap >= -tolerance*max(1,abs(norms(2))), ...
            'TSF:Precision','Negative pooled-neutrality difference.');
        pooledNeutrality = exp(-norms(3))*(-expm1(-max(gap,0)));
        membershipLoss = exp(-norms(1));
        pooledMembership = -expm1(-norms(1));
        pooledNonmembership = exp(-norms(2));

        logMembership = log1p(-membershipLoss);
        logCumulative = log1p(pooledNeutrality-membershipLoss);
        assert(logCumulative <= tolerance, ...
            'TSF:Closure','Invalid pooled cumulative component.');
        logCumulative = min(logCumulative,0);
        sAggregate = exp(logMembership/scale);
        dAggregate = -expm1(log1p(-pooledNonmembership)/scale);
        cumulativeAggregate = exp(logCumulative/scale);
        if pooledNeutrality == 0
            hAggregate = 0;
        elseif pooledMembership > 0
            hAggregate = cumulativeAggregate*(-expm1( ...
                -log1p(pooledNeutrality/pooledMembership)/scale));
        else
            hAggregate = cumulativeAggregate;
        end
        aggregatePowers = [sAggregate,dAggregate,hAggregate];
        assert(all(isfinite(aggregatePowers)) && all(aggregatePowers >= -tolerance) && ...
            sum(aggregatePowers) <= 1+tolerance, ...
            'TSF:Precision','Aggregate failed the numerical feasibility check.');

        % 4.7 转回 (S,D,H)，计算评分、精确度和确定度。
        aggregate(a,:) = aggregatePowers.^(1/gamma);
        indices(a,:) = [(1+sAggregate-dAggregate)/2, ...
                       sAggregate+dAggregate, sAggregate+dAggregate+hAggregate];
    end

    % 4.8 依次按评分、精确度、确定度降序排列。
    [~,ranking] = sortrows(-indices,[1,2,3]);
    for j = 1:nAlternatives-1
        a = ranking(j);
        b = ranking(j+1);
        if abs(indices(a,1)-indices(b,1)) <= 1e-12
            assert(isequal(mass(a,:,1:2),mass(b,:,1:2)), ...
                'TSF:Precision','Unresolved near-primary comparison.');
            if abs(indices(a,3)-indices(b,3)) <= 1e-12
                assert(isequal(mass(a,:,:),mass(b,:,:)), ...
                    'TSF:Precision','Unresolved near-neutrality comparison.');
            end
        end
    end

    allAggregates(:,:,setting) = aggregate;
    allIndices(:,:,setting) = indices;
    caseTables{setting} = table(alternatives,aggregate(:,1),aggregate(:,2),aggregate(:,3), ...
        indices(:,1),indices(:,2),indices(:,3), ...
        'VariableNames',{'Alternative','S','D','H','Score','Accuracy','Determination'});
    caseRankings{setting} = strjoin(cellstr(alternatives(ranking)),' > ');

    fprintf('\n%s weights\n',settingNames{setting});
    fprintf('Alternative      S        D        H      Score   Accuracy   Determination\n');
    for a = 1:nAlternatives
        fprintf('%-10s   %.4f   %.4f   %.4f   %.4f    %.4f       %.4f\n', ...
            char(alternatives(a)),aggregate(a,:),indices(a,:));
    end
    fprintf('Ranking: %s\n',caseRankings{setting});
end

%% 5. 保留工作区结果，并导出两组 CSV
equalAggregate = allAggregates(:,:,1);
equalIndices = allIndices(:,:,1);
unequalAggregate = allAggregates(:,:,2);
unequalIndices = allIndices(:,:,2);
scriptDir = fileparts(mfilename('fullpath'));
outputDir = fullfile(scriptDir,'results', ...
    ['case_' char(datetime('now','Format','yyyyMMdd_HHmmss_SSS'))]);
assert(~isfolder(outputDir),'TSF:Output','Output directory already exists.');
mkdir(outputDir);
for setting = 1:nSettings
    writetable(caseTables{setting},fullfile(outputDir, ...
        [settingNames{setting} '_weights.csv']),'Encoding','UTF-8');
end
fprintf('\nComputed tables saved to %s\n',outputDir);
