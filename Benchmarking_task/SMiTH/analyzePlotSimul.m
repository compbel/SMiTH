clear;
nTests = 500;
minSize = 5;
maxSize = 30;

load('nVertTest500.mat');

netPriors = {'True','SF','Prufer'};
netPriorsNames = {'Degenerate','Scale-free','Uniform'};
percSampDP = [75 99 0];
plotData = struct('groupsTests',{},'groupsSamples',{},'groupsTestsExt1',{},'groupsPercentiles',{},...
    'percHomomorAll',{},'fscoreAll',{},'fscoreSpanAll',{},'percentileTrue',{});


models = {'noObj_noConvex_noSampPars','obj_convex_noSampParsDP','obj_convex_sampParsObjDP',... 
    'obj_convex_sampParsDP'};
model_names = {'Unconstrained','Convex','Convex and max compact','Convex and compact'};

for net = 1:length(netPriors)
    percHomomorAll = [];
    fscoreAll = [];
    fscoreSpanAll = [];
    groupsTests = [];
    groupsTestsAll = [];
    groupsSamples = [];
    groupsPercentiles = [];
    histAll = [];
    histAllTrue = [];
    percentileTrue = [];
    precisionThr = [];
    fOptThr = [];
    precisionThrW = [];
    fOptThrW = [];
    aucAll = [];
    numFeas = -ones(length(models),nTests);
    numSamp = -ones(length(models),nTests);
    
    fscoreSpanAllArray = -ones(length(models),nTests);
    
    for i =  1:length(models)
        if ~(startsWith(models{i},'obj') | startsWith(models{i},'noObj'))
            continue;
        end
        resFolder = ['results' filesep models{i} filesep netPriors{net}];
        percHomomor = -ones(1,nTests);
        fscoreTest = [];
        sizes = [];
        for t = 1:nTests
            t
            filename = [resFolder filesep 'res_t' num2str(t) '.mat'];
            if ~exist([resFolder filesep 'res_t' num2str(t) '.mat'],'file')
                continue;
            end
            if ~ismember(nVertTest(t),minSize:maxSize) %|| (numViol(t) > 0)
                continue;
            end
            load(filename);
            n = size(AMTNtrue,1);
            nSamp = length(sensSamp);
            numSamp(i,t) = nSamp;
            ind = (sensSamp ~= -1);
            percHomomor(t) = sum(ind)/nSamp;
            numFeas(i,t) = sum(ind)/nSamp;
            fscoreTest = [fscoreTest sensSamp(ind)];
            AMTNinferSamp = AMTNinferSamp(ind);
   
    
            AMCons = sum(cat(3, AMTNinferSamp{:}), 3)/sum(ind);
            edges = 0:0.1:1;
            upperTriIndex = triu(true(size(AMCons)),1);
    
            if isempty(AMCons)
                continue;
            end
    
    
            AMspan = full(adjacency(minspantree(graph(-AMCons))));
            fscoreSpan = sum(sum(AMspan.*AMTNtrue))/sum(sum(AMTNtrue));
            fscoreSpanAllArray(i,t) = fscoreSpan;
            fscoreSpanAll = [fscoreSpanAll fscoreSpan];
            comparMatr = (AMCons(upperTriIndex) <= (AMCons(AMTNtrue==1))');
            percentileTrue = [percentileTrue sum(comparMatr,1)/numel(AMCons(upperTriIndex))];
            groupsPercentiles = [groupsPercentiles i*ones(1,length(sum(comparMatr,1)/numel(AMCons(upperTriIndex))))];
            thr = min(AMCons(AMTNtrue==1));
            AMConsThr = AMCons.*(AMCons>=thr);
            precisionThr = [precisionThr sum(sum(AMConsThr.*AMTNtrue))/sum(sum(AMConsThr))];
            precisionThrW = [precisionThrW sum(sum(AMCons.*AMConsThr.*AMTNtrue))/sum(sum(AMCons.*AMConsThr))];
            sizes = [sizes length(AMspan)];
            thr = 0:0.01:1;
            best_f = 0;
            best_f_w = 0;
            tpr = zeros(1,length(thr));
            fpr = zeros(1,length(thr)); 
            prec = zeros(1,length(thr));
            rec = zeros(1,length(thr)); 
            for th = 1:length(thr)
                AMTNinfer = (AMCons >= thr(th));
                fscore = sum(sum(AMTNinfer.*AMTNtrue))/sum(sum(AMTNtrue));
                if sum(sum(AMTNinfer)) > 0
                    precision = sum(sum(AMTNinfer.*AMTNtrue))/(sum(sum(AMTNinfer))/2);
                else
                    precision = 0;
                end
                f = 2*fscore*precision/(fscore+precision);
                if f > best_f
                    best_f = f;
                end
                fscore_w = sum(sum(AMCons.*AMTNinfer.*AMTNtrue))/sum(sum(AMCons.*AMTNtrue));
                precision_w = sum(sum(AMCons.*AMTNinfer.*AMTNtrue))/sum(sum(AMCons.*AMTNinfer));
                f_w = 2*fscore_w*precision_w/(fscore_w+precision_w);
                if f_w > best_f_w
                    best_f_w = f_w;
                end
                tpr(th) = fscore;
                AMTNfalse = 1 - (AMTNtrue+AMTNtrue'>0) - eye(n,n);
                fpr(th) = sum(sum(AMTNinfer.*AMTNfalse))/(n*(n-1)-2*sum(sum(AMTNtrue)));
                rec(th) = fscore;
                prec(th) = precision;
            end
            fOptThr = [fOptThr best_f];
            fOptThrW = [fOptThrW best_f_w];
            auc = trapz(flip([rec 0]'),flip([prec 1]'));
            aucAll = [aucAll auc];
            aucAllArray(i,t) = auc;


        end
        ind = find(percHomomor > -1);
        percHomomorAll = [percHomomorAll percHomomor(ind)];
        groupsTestsAll = [groupsTestsAll i*ones(1,length(ind))];

        ind = find(percHomomor > 0);
        fscoreAll = [fscoreAll fscoreTest];
        groupsTests = [groupsTests i*ones(1,length(ind))];
        groupsSamples = [groupsSamples i*ones(1,length(fscoreTest))];
    end
    
    groupsTestsExt = groupsTests;
    i = find(strcmp(models, 'weighted'));
    if ~isempty(i)
        groupsTestsExt = groupsTests;
        resFolder = ['results' filesep models{i} filesep netPriors{net}];
        for t = 1:nTests
            t
            filename = [resFolder filesep 'res_t' num2str(t) '.mat'];
            if ~exist([resFolder filesep 'res_t' num2str(t) '.mat'],'file')
                continue;
            end
            load(filename);
            fscoreSpanAll = [fscoreSpanAll fscore];
            groupsTestsExt = [groupsTestsExt i];
        end
    end
     
    groupsTestsExt1 = groupsTestsExt;
    
    i = find(strcmp(models, 'TNet'));
    if ~isempty(i)
        tnet_file = ['Results for TNet 500' filesep 'tnet_500_MTS_results.csv'];
        tnet_res = readmatrix(tnet_file);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (tnet_res(ind,2))'];
        fscoreSpanAllArray(i,:) = (tnet_res(:,2))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end
    
    i = find(strcmp(models, 'Cassiopeia'));
    if ~isempty(i)
        cassiop_file = ['Cassiopeia_all_results' filesep 'cassiopeia_fscore_for_MST_results.csv'];
        cassiop_res = readmatrix(cassiop_file);
        cassiop_res = cassiop_res(1:nTests,:);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (cassiop_res(ind,2))'];
        fscoreSpanAllArray(i,:) = (cassiop_res(:,2))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end
    
    
    % % PS (parallel single-source seeding)	Each metastatic site is seeded from only one other anatomical site, i.e. G is a multi-tree.
    % % M (multi-source seeding)	A metastatic site may be seeded from multiple anatomical sites, but no directed cycles are introduced. That is, G is multi-DAG.
    % %R (reseeding)	Directed cycles in G are allowed.
    % % S (single-source seeding) tree
    
    i = find(strcmp(models, 'Machina'));
    if ~isempty(i)
        machina_file = ['Machina_results' filesep 'machina_fscoreitivity_R.csv'];
        machina_res = readmatrix(machina_file);
        machina_res = machina_res(1:nTests,:);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (machina_res(ind,2))'];
        fscoreSpanAllArray(i,:) = (machina_res(:,2))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end
    
    i = find(strcmp(models, 'Phyloscanner'));
    if ~isempty(i)
        % phylos_file = ['Phyloscanner_results' filesep 'all_stats_phyloscanner_new_k0.csv'];
        phylos_file = ['Phyloscanner_results' filesep 'all_stats_phyloscanner_new_adjacency_goes_through_undefined_k0.csv'];
        phylos_res = readmatrix(phylos_file);
        phylos_res = phylos_res(1:nTests,:);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (phylos_res(ind,4))'];
        fscoreSpanAllArray(i,:) = (phylos_res(:,4))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end

    plotData(net).groupsTests = groupsTests;
    plotData(net).groupsSamples = groupsSamples;
    plotData(net).groupsTestsExt1 = groupsTestsExt1;
    plotData(net).groupsPercentiles = groupsPercentiles;
    plotData(net).percHomomorAll = percHomomorAll;
    plotData(net).fscoreAll = fscoreAll;
    plotData(net).fscoreSpanAll = fscoreSpanAll;
    plotData(net).percentileTrue = percentileTrue;
    plotData(net).numSamp = numSamp;
    plotData(net).aucAll = aucAll;
    plotData(net).groupsTestsAll = groupsTestsAll;
end

%% 
% percent compatible
nRow = 1;
nCol = length(netPriors);
figure
[ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 1800, 570])
% set(gcf, 'Position',  [100, 100, 2000, 870]) %for RECOMB
for net = 1:length(netPriors)
    axes(ha(net));
    i = find(strcmp(model_names,'Convex and max compact'));
    ind = (plotData(net).groupsTestsAll ~= i);
    groupsTests_curr = plotData(net).groupsTestsAll(ind);
    groupsTests_curr(groupsTests_curr > i) = groupsTests_curr(groupsTests_curr > i) - 1;
    bc = boxchart(groupsTests_curr,plotData(net).percHomomorAll(ind),'GroupByColor',groupsTests_curr,...
        "BoxWidth",3,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    bc(3).SeriesIndex = 4;
    ylabel('percent','FontSize', 18,'FontWeight','bold');
    ax = gca;
    ax.XLim = [0 4];
    ax.FontSize = 16;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
    text(0.99*xl(2),0.99*yl(2),['(' char(97+net-1) ')'],'FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
    title(netPriorsNames{net},'FontSize', 18);
end
lgd = legend(model_names(setdiff(1:length(model_names),i)),'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.36;
lgd.Position(2) = 0.03;

%% 
% f-score sample
nRow = 1;
nCol = length(netPriors);
figure
[ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 1800, 570])
for net = 1:length(netPriors)
    axes(ha(net));
    i = find(strcmp(model_names,'Convex and max compact'));
    ind = (plotData(net).groupsSamples ~= i);
    groupsSamples_curr = plotData(net).groupsSamples(ind);
    groupsSamples_curr(groupsSamples_curr > i) = groupsSamples_curr(groupsSamples_curr > i) - 1;
    bc = boxchart(groupsSamples_curr,plotData(net).fscoreAll(ind),'GroupByColor',groupsSamples_curr,...
        "BoxWidth",3,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    bc(3).SeriesIndex = 4;
    ylabel('f-score','FontSize', 18,'FontWeight','bold');
    ax = gca;
    ax.XLim = [0 4];
    ax.FontSize = 16;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
    text(0.99*xl(2),0.99*yl(2),['(' char(97+net-1) ')'],'FontSize',40,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
    title(netPriorsNames{net},'FontSize', 18);
end
lgd = legend(model_names(setdiff(1:length(model_names),i)),'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.36;
lgd.Position(2) = 0.03;
exportgraphics(gcf,['figures' filesep 'f_sample.png'],'Resolution',600)
% exportgraphics(gcf,['figures' filesep 'f_sample1.png'],'Resolution',600) %for RECOMB

%% 
%plot percentiles
nRow = 1;
nCol = length(netPriors);
figure
[ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 1800, 570])
for net = 1:length(netPriors)
    axes(ha(net));
    i = find(strcmp(model_names,'Convex and max compact'));
    ind = (plotData(net).groupsPercentiles ~= i);
    groupsPercentiles_curr = plotData(net).groupsPercentiles(ind);
    groupsPercentiles_curr(groupsPercentiles_curr > i) = groupsPercentiles_curr(groupsPercentiles_curr > i) - 1;
    bc = boxchart(groupsPercentiles_curr,plotData(net).percentileTrue(ind),'GroupByColor',groupsPercentiles_curr,...
        "BoxWidth",3,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    bc(3).SeriesIndex = 4;
    ylabel('f-score','FontSize', 18,'FontWeight','bold');
    title(netPriorsNames{net},'FontSize', 18);
    ax = gca;
    ax.XLim = [0 4];
    ax.FontSize = 16;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
    text(0.99*xl(2),0.1*yl(2),['(' char(100+net-1) ')'],'FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')    
end
lgd = legend(model_names(setdiff(1:length(model_names),i)),'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.35;
lgd.Position(2) = 0.03;
exportgraphics(gcf,['figures' filesep 'percentile_true.png'],'Resolution',600)
%% 

%plot auc
nRow = 1;
nCol = length(netPriors);
figure
[ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 1800, 570])
% set(gcf, 'Position',  [100, 100, 2000, 870]) %for RECOMB
for net = 1:length(netPriors)
    axes(ha(net));
    i = find(strcmp(model_names,'Convex and max compact'));
    ind = (plotData(net).groupsTestsExt1 ~= i);
    groupsAucAll_curr = plotData(net).groupsTestsExt1(ind);
    groupsAucAll_curr(groupsAucAll_curr > i) = groupsAucAll_curr(groupsAucAll_curr > i) - 1;
    bc = boxchart(groupsAucAll_curr,plotData(net).aucAll(ind),'GroupByColor',groupsAucAll_curr,...
        "BoxWidth",3,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    bc(3).SeriesIndex = 4;
    ylabel('AUC','FontSize', 18,'FontWeight','bold');
    ax = gca;
    ax.FontSize = 16;
    % ax.FontSize = 30; %for RECOMB
    ax.FontWeight = 'bold';
    ax.XLim = [0 4];
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
    text(0.99*xl(2),0.1*yl(2),['(' char(100+net-1) ')'],'FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
    title(netPriorsNames{net},'FontSize', 18);
end
lgd = legend(model_names(setdiff(1:length(model_names),i)),'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.35;
lgd.Position(2) = 0.03;
exportgraphics(gcf,['figures' filesep 'auc.png'],'Resolution',600)

%% 

% spanning tree

nRow = 1;
nCol = length(netPriors);
figure
[ha, pos] = tight_subplot(nRow,nCol,[.045 .03],[.1 0.05],[.1 .1]);
set(gcf, 'Position',  [100, 100, 1800, 570])
for net = 1:length(netPriors)
    axes(ha(net));
    i = find(strcmp(model_names,'Convex and max compact'));
    ind = (plotData(net).groupsTestsExt1 ~= i);
    groupsTestsExt1_curr = plotData(net).groupsTestsExt1(ind);
    groupsTestsExt1_curr(groupsTestsExt1_curr > i) = groupsTestsExt1_curr(groupsTestsExt1_curr > i) - 1;
    bc = boxchart(groupsTestsExt1_curr,plotData(net).fscoreSpanAll(ind),'GroupByColor',groupsTestsExt1_curr,...
        "BoxWidth",3,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    bc(3).SeriesIndex = 4;
    ylabel('f-score','FontSize', 18,'FontWeight','bold');
    title(netPriorsNames{net},'FontSize', 18);
    ax = gca;
    ax.XLim = [0 4];
    ax.FontSize = 16;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
    text(0.99*xl(2),0.1*yl(2),['(' char(97+net-1) ')'],'FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
end
lgd = legend(model_names(setdiff(1:length(model_names),i)),'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.35;
lgd.Position(2) = 0.03;
exportgraphics(gcf,['figures' filesep 'fscore_span.png'],'Resolution',600)
%% 
% algorithm comparision
clear;
nTests = 500;
minSize = 5;
maxSize = 30;

load('nVertTest500.mat');

netPriors = {'SF','Prufer'};
netPriorsNames = {'Scale-free','Uniform'};
percSampDP = [75 99 0];
plotData = struct('groupsTestsExt1',{},'fscoreSpanAll',{});


models = {'obj_convex_sampParsObjDP','Cassiopeia','Stratus','MACHINA','Phyloscanner','TNet'};
model_names = {'SMITH SF','SMITH Uniform','Cassiopeia','STraTUS','MACHINA','Phyloscanner','TNet'};

fscoreSpanAll = [];
groupsTests = [];
fscoreSpanAllArray = -ones(length(models),nTests);
runTimeAll = cell(length(netPriors),nTests);
j = 1;
for net = 1:length(netPriors)
    for i =  1:length(models)
        if ~(startsWith(models{i},'obj') | startsWith(models{i},'noObj'))
            continue;
        end
        resFolder = ['results' filesep models{i} filesep netPriors{net}];
        percHomomor = -ones(1,nTests);
        for t = 1:nTests
            t
            filename = [resFolder filesep 'res_t' num2str(t) '.mat'];
            if ~exist([resFolder filesep 'res_t' num2str(t) '.mat'],'file')
                continue;
            end
            if ~ismember(nVertTest(t),minSize:maxSize) %|| (numViol(t) > 0)
                continue;
            end
            load(filename);
            n = size(AMTNtrue,1);
            nSamp = length(sensSamp);
            numSamp(i,t) = nSamp;
            ind = (sensSamp ~= -1);
            percHomomor(t) = sum(ind)/nSamp;
            AMTNinferSamp = AMTNinferSamp(ind);
    
            if strcmp(models{i}, 'obj_convex_sampParsObjDP')
                objSamp = objSamp(ind);
                p = prctile(objSamp, percSampDP(net));
                ind1 = (objSamp >= p);
                AMTNinferSamp = AMTNinferSamp(ind1);
            end
    
            AMCons = sum(cat(3, AMTNinferSamp{:}), 3)/sum(ind);
            upperTriIndex = triu(true(size(AMCons)),1);

            if isempty(AMCons)
                continue;
            end
    
    
            AMspan = full(adjacency(minspantree(graph(-AMCons))));
            fscoreSpan = sum(sum(AMspan.*AMTNtrue))/sum(sum(AMTNtrue));
            fscoreSpanAllArray(j,t) = fscoreSpan;
            fscoreSpanAll = [fscoreSpanAll fscoreSpan];
            runTimeAll{net,t} = timeSamp(timeSamp ~= -1);
        end
        ind = find(percHomomor > 0);
        percHomomor = percHomomor(ind);
        groupsTests = [groupsTests j*ones(1,length(ind))];
        j = j + 1;
    end
end
    
    groupsTestsExt1 = groupsTests;
         
    i = find(strcmp(model_names, 'TNet'));
    if ~isempty(i)
        tnet_file = ['Results for TNet 500' filesep 'tnet_500_MTS_results.csv'];
        tnet_res = readmatrix(tnet_file);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (tnet_res(ind,2))'];
        fscoreSpanAllArray(i,:) = (tnet_res(:,2))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end
    
    i = find(strcmp(model_names, 'Cassiopeia'));
    if ~isempty(i)
        cassiop_file = ['Cassiopeia_all_results' filesep 'cassiopeia_sens_for_MST_results.csv'];
        cassiop_res = readmatrix(cassiop_file);
        cassiop_res = cassiop_res(1:nTests,:);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (cassiop_res(ind,2))'];
        fscoreSpanAllArray(i,:) = (cassiop_res(:,2))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end
    
    
    % % PS (parallel single-source seeding)	Each metastatic site is seeded from only one other anatomical site, i.e. G is a multi-tree.
    % % M (multi-source seeding)	A metastatic site may be seeded from multiple anatomical sites, but no directed cycles are introduced. That is, G is multi-DAG.
    % %R (reseeding)	Directed cycles in G are allowed.
    % % S (single-source seeding) tree

    
    i = find(strcmp(model_names, 'MACHINA'));
    if ~isempty(i)
        machina_file = ['Machina_results' filesep 'all_stats_machina_new_S.csv'];
        machina_res = readmatrix(machina_file);
        machina_res = machina_res(1:nTests,:);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (machina_res(ind,2))'];
        fscoreSpanAllArray(i,:) = (machina_res(:,2))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end
    
    i = find(strcmp(model_names, 'Phyloscanner'));
    if ~isempty(i)
        phylos_file = ['Phyloscanner_results' filesep 'all_stats_phyloscanner_new_adjacency_goes_through_undefined_k0.csv'];
        phylos_res = readmatrix(phylos_file);
        phylos_res = phylos_res(1:nTests,:);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (phylos_res(ind,4))'];
        fscoreSpanAllArray(i,:) = (phylos_res(:,4))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end

    i = find(strcmp(model_names, 'STraTUS'));
    if ~isempty(i)
        alg_file = ['Stratus_results' filesep 'all_stats_stratus_1000_samples.csv'];
        alg_res = readmatrix(alg_file);
        alg_res = alg_res(1:nTests,:);
        ind = (nVertTest>=5)&(nVertTest<=maxSize);
        fscoreSpanAll = [fscoreSpanAll (alg_res(ind,4))'];
        fscoreSpanAllArray(i,:) = (alg_res(:,4))';
        groupsTestsExt1 = [groupsTestsExt1 i*ones(1,sum(ind))];
    end

%% 
% spanning tree
figure
set(gcf, 'Position',  [100, 100, 1300, 570])
    % axes(ha(net));
    bc = boxchart(groupsTestsExt1,fscoreSpanAll,'GroupByColor',groupsTestsExt1,...
        "BoxWidth",6,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    ylabel('f-score','FontSize', 18,'FontWeight','bold');
    ax = gca;
    ax.XLim = [0 8];
    ax.FontSize = 16;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
    text(0.99*xl(2),0.99*yl(2),'(a)','FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
lgd = legend(model_names,'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.14;
lgd.Position(2) = 0.03;
exportgraphics(gcf,['figures' filesep 'fscore_span.png'],'Resolution',600)
%% 
% accuracy and running time by size
sizeBoxes = minSize:3:maxSize;
sizeBoxes(end) = sizeBoxes(end)+1;
fscoreSpanSize = zeros(size(fscoreSpanAllArray,1),length(sizeBoxes)-1);
numSize = zeros(1,length(sizeBoxes)-1);
for i = 1:size(fscoreSpanAllArray,1)
    for j = 1:length(sizeBoxes)-1
        ind = find(nVertTest>=sizeBoxes(j) & nVertTest<sizeBoxes(j+1));
        numSize(j) = length(ind);
        fscoreSpanSize(i,j) = median(fscoreSpanAllArray(i,ind));
    end
end

runtimeSpanSize = zeros(size(runTimeAll,1),length(sizeBoxes)-1);
for i = 1:size(runTimeAll,1)
   for j = 1:length(sizeBoxes)-1
        ind = find(nVertTest>=sizeBoxes(j) & nVertTest<sizeBoxes(j+1));
        runtimeSpanSize(i,j) = median(horzcat(runTimeAll{i,ind}));
    end
end

xsizes = sizeBoxes;
xsizes(end) = 30;
figure
set(gcf, 'Position',  [100, 100, 500, 570])
plot(sizeBoxes(2:end),fscoreSpanSize','LineWidth',4);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
xl = xlim;
yl = ylim;
text(0.99*xl(2),1.1*yl(2),'(b)','FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
lgd = legend(model_names,'Location','southoutside','FontSize',12,'FontWeight','bold');
lgd.NumColumns = round(length(model_names)/2);
xlim([sizeBoxes(2) sizeBoxes(end)]);
ylim([0 1])
xlabel('number of populations','FontSize', 18,'FontWeight','bold');
ylabel('f-score','FontSize', 18,'FontWeight','bold');
exportgraphics(gcf,['figures' filesep 'fscore_time.png'],'Resolution',600)

figure
set(gcf, 'Position',  [100, 100, 500, 570])
plot(sizeBoxes(2:end),runtimeSpanSize','LineWidth',4);
ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
xl = xlim;
yl = ylim;
text(0.99*xl(2),1.15*yl(2),'(c)','FontSize',30,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','top')
lgd = legend(netPriorsNames,'Location','southoutside','FontSize',12,'FontWeight','bold');
lgd.NumColumns = length(netPriorsNames);
xlim([sizeBoxes(2) sizeBoxes(end)]);
ylim([0 0.8])
xlabel('number of populations','FontSize', 18,'FontWeight','bold');
ylabel('time per candidate tree (sec)','FontSize', 18,'FontWeight','bold');
exportgraphics(gcf,['figures' filesep 'runtime_time.png'],'Resolution',600)

%% 
% varied sample sizes
clear;
numSeq = [100 50 10 1];
model_names = {'S = 100', 'S = 50', 'S = 10', 'S = 1'};
groups = [];
fscores = [];
for s = 1:length(numSeq)
    s = 3;
    resFolder = ['results' filesep 'obj_convex_sampParsObjDP' filesep 'SF_exp' int2str(numSeq(s))];
    nTests = 500;
    nVertTestRev = zeros(1,nTests);
    sampSize = zeros(1,nTests);
    rec = zeros(1,nTests);
            for t = 1:nTests
                t
                filename = [resFolder filesep 'res_t' num2str(t) '.mat'];
                if ~exist(filename,'file')
                    continue;
                end
                load(filename);
                nVertTestRev(t) = length(AMTNtrue);
                sampSize(t) = sum(sensSamp ~= -1);
                ind = (sensSamp ~= -1);
                AMTNinferSamp = AMTNinferSamp(ind);
                objSamp = objSamp(ind);
                p = prctile(objSamp, 75);
    
                ind1 = (objSamp >= p);
    
                AMTNinferSamp = AMTNinferSamp(ind1);
                AMCons = sum(cat(3, AMTNinferSamp{:}), 3)/sum(ind);
                upperTriIndex = triu(true(size(AMCons)),1);
                AMspan = full(adjacency(minspantree(graph(-AMCons))));
                recallSpan = sum(sum(AMspan.*AMTNtrue))/sum(sum(AMTNtrue));
                rec(t) = recallSpan;
            end
    
            ind = nVertTestRev >= 5 & nVertTestRev<=30;
            mean(rec(ind))
            median(rec(ind))
            std(rec(ind))
            groups = [groups s*ones(1,sum(ind))];
            fscores = [fscores rec(ind)];
            sum(ind)
end

figure
set(gcf, 'Position',  [100, 100, 1300, 570])
    bc = boxchart(groups,fscores,'GroupByColor',groups,...
        "BoxWidth",3,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    ylabel('f-score','FontSize', 18,'FontWeight','bold');
    % title('Algorithms comparision','FontSize', 18);
    ax = gca;
    ax.XLim = [0 5];
    ax.FontSize = 16;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
lgd = legend(model_names,'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.35;
lgd.Position(2) = 0.03;
exportgraphics(gcf,['figures' filesep 'samp_size.png'],'Resolution',600)

%% 
% random transmission rates
clear;
resFolder = ['results' filesep 'obj_convex_sampParsObjDP' filesep 'SF_randTransRates'];
nTests = 500;
nVertTestRev = zeros(1,nTests);
sampSize = zeros(1,nTests);
rec = zeros(1,nTests);
        for t = 1:nTests
            t
            filename = [resFolder filesep 'res_t' num2str(t) '.mat'];
            if ~exist(filename,'file')
                continue;
            end
            load(filename);
            nVertTestRev(t) = length(AMTNtrue);
            sampSize(t) = sum(sensSamp ~= -1);
            ind = (sensSamp ~= -1);
            AMTNinferSamp = AMTNinferSamp(ind);
            objSamp = objSamp(ind);
            p = prctile(objSamp, 75);

            ind1 = (objSamp >= p);

            AMTNinferSamp = AMTNinferSamp(ind1);
            AMCons = sum(cat(3, AMTNinferSamp{:}), 3)/sum(ind);
            upperTriIndex = triu(true(size(AMCons)),1);
            AMspan = full(adjacency(minspantree(graph(-AMCons))));
            recallSpan = sum(sum(AMspan.*AMTNtrue))/sum(sum(AMTNtrue));
            rec(t) = recallSpan;
        end

        ind = nVertTestRev >= 5 & nVertTestRev<=30;
        sum(ind)
        mean(rec(ind))
        median(rec(ind))
        std(rec(ind))
%% 
% random sample sizes
clear;
resFolder = ['results' filesep 'obj_convex_sampParsObjDP' filesep 'SF_undersamp'];
nTests = 500;
nVertTestRev = zeros(1,nTests);
sampSize = zeros(1,nTests);
sampPat = zeros(1,nTests);
rec = zeros(1,nTests);
prec = zeros(1,nTests);
        for t = 1:nTests
            t
            filename = [resFolder filesep 'res_t' num2str(t) '.mat'];
            if ~exist(filename,'file')
                continue;
            end
            load(filename);
            nVertTestRev(t) = length(AMTNtrue);
            sampSize(t) = sum(sensSamp ~= -1);
            sampPat(t) = retainPat;
            ind = (sensSamp ~= -1);
            AMTNinferSamp = AMTNinferSamp(ind);
            objSamp = objSamp(ind);
            p = prctile(objSamp, 75);

            ind1 = (objSamp >= p);

            AMTNinferSamp = AMTNinferSamp(ind1);
            AMCons = sum(cat(3, AMTNinferSamp{:}), 3)/sum(ind);
            upperTriIndex = triu(true(size(AMCons)),1);
            AMspan = full(adjacency(minspantree(graph(-AMCons))));
            recallSpan = sum(sum(AMspan.*AMTNtrue))/sum(sum(AMTNtrue));
            precSpan = sum(sum(AMspan.*AMTNtrue))/(sum(sum(AMspan))/2);
            rec(t) = recallSpan;
            prec(t) = precSpan;
        end

        fscore = 2*(prec.*rec)./(prec+rec);

        ind = nVertTestRev >= 5 & nVertTestRev<=30;
        mean(fscore(ind))
        median(fscore(ind))
        std(fscore(ind))
        
%% 
% random patients undersampling
clear;
resFolder = ['results' filesep 'obj_convex_sampParsObjDP' filesep 'SF_undersampPat'];
nTests = 500;
nVertTestRev = zeros(1,nTests);
sampSize = zeros(1,nTests);
sampPat = zeros(1,nTests);
rec = zeros(1,nTests);
prec = zeros(1,nTests);
        for t = 1:nTests
            t
            filename = [resFolder filesep 'res_t' num2str(t) '.mat'];
            if ~exist(filename,'file')
                continue;
            end
            load(filename);
            nVertTestRev(t) = length(AMTNtrue);
            sampSize(t) = sum(sensSamp ~= -1);
            sampPat(t) = retainPat;
            ind = (sensSamp ~= -1);
            AMTNinferSamp = AMTNinferSamp(ind);
            objSamp = objSamp(ind);
            p = prctile(objSamp, 75);

            ind1 = (objSamp >= p);

            AMTNinferSamp = AMTNinferSamp(ind1);
            AMCons = sum(cat(3, AMTNinferSamp{:}), 3)/sum(ind);
            upperTriIndex = triu(true(size(AMCons)),1);
            AMspan = full(adjacency(minspantree(graph(-AMCons))));
            recallSpan = sum(sum(AMspan.*AMTNtrue))/sum(sum(AMTNtrue));
            precSpan = sum(sum(AMspan.*AMTNtrue))/(sum(sum(AMspan))/2);
            rec(t) = recallSpan;
            prec(t) = precSpan;
        end

        prec(isnan(prec)) = 0;
        rec(isnan(rec)) = 0;
        fscore = 2*(prec.*rec)./(prec+rec);
        fscore(isnan(fscore)) = 0;
        ind = nVertTestRev >= 15 & nVertTestRev<=30;
        mean(fscore(ind))
        median(fscore(ind))
        std(fscore(ind))
%% 

% plot random sampling and random transmission rates
clear;
fscores = [];
groups = [];
model_names = {'basic', 'random migration rates', 'random sampling rates'};

load("results/fscores_basic.mat");
fscores = [fscores rec];
groups = [groups ones(1,length(rec))];
load("results/fscores_varrate.mat");
fscores = [fscores rec];
groups = [groups 2*ones(1,length(rec))];
load("results/fscores_varsamp.mat");
fscores = [fscores fscore];
groups = [groups 3*ones(1,length(fscore))];


figure
set(gcf, 'Position',  [100, 100, 1300, 570])
    % axes(ha(net));
    bc = boxchart(groups,fscores,'GroupByColor',groups,...
        "BoxWidth",3,'BoxFaceAlpha',0.5,'BoxMedianLineColor','k');
    ylabel('f-score','FontSize', 18,'FontWeight','bold');
    ax = gca;
    ax.XLim = [0 4];
    ax.FontSize = 16;
    ax.FontWeight = 'bold';
    set(gca,'XTickLabel',[]);
    ylim([0 1]);
    grid on
    xl = xlim;
    yl = ylim;
lgd = legend(model_names,'Location','southoutside','FontSize',18,'FontWeight','bold');
lgd.NumColumns = length(model_names);
lgd.Position(1) = 0.35;
lgd.Position(2) = 0.03;
exportgraphics(gcf,['figures' filesep 'misc_rand.png'],'Resolution',600)