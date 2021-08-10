function SCHMM_plot_results(chr,maxCN,datafile,resultsfile,plotsdir,barcode)
%this function is used to plot aberration detection results of a chromosome
%----------------------read results  ------------------------%

cns = 0.01;
mcns = 0.5;
for i = 1:maxCN
	k = max(floor(i/2), floor((i+1)/2));
	for j = k:i
		cns = [cns i];
		mcns = [mcns j];
	end
end

fid = fopen(resultsfile,'r');
if fid == -1
    error(['Can not open result file: ' resultsfile]);
end

%get estimated global parameters from the first row of the result file
lambda = [];
mu = [];

while 1
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'StartPos')),break,end
    %Lambda
    result1 = regexp(tline,'Copy neutral read counts:\s*(\S+)','tokens','once');
    if ~isempty(result1)
        lambda = str2double(result1{1});
    end
    %BAF mu
    if ~isempty(strfind(tline,'Mus of BBDs for major allele depth'))
        tline = fgetl(fid);
        results = regexp(strtrim(tline), '\s', 'split');
        mu = str2double(results);
    end
end
%report errors if these values are not parsed successfully
if isempty(lambda)
    error(['Cannot read estimated copy neutral read counts from ',resultsfile]);
end
if isempty(mu)
    error(['Cannot read estimated major-allele frequency from ',resultsfile]);
end

%then read the results
results = textscan(fid, '%f%f%f%f%f%*s%f', 'treatAsEmpty', {'NA', 'na'});
fclose(fid);
chr_seg = results{1};
pstart_seg = results{2};
pend_seg = results{3};
cn_seg = results{4};
mcn_seg = results{5};
% AI_seg = results{6};
score_seg = results{6};
clear results;

%load data
fid = fopen(datafile,'r');
if fid == -1
    error(['Can not open data file: ' datafile]);
end
results = textscan(fid, repmat('%f', 1, 7), 'HeaderLines', 1);
fclose(fid);
data_chr_all = results{1};
data_pos_all = results{2};
data_bd_all = results{3};
data_td_all = results{4};
data_rc_all = results{5};
clear results;

%--------------- plot figures ---------------------%
rc_colors = [0.5 0.5 0.5;
             0 0.9 0;
             0 0 0.9;
             0.9 0 0];
baf_colors = [0 250 0;
             0 0 250;
             250 0 0]./255;
         
if chr ~= 0
    tv = ismember(data_chr_all,chr);
    data_rc = data_rc_all(tv);
    data_bd = data_bd_all(tv);
    data_td = data_td_all(tv);
    data_baf = data_bd./data_td;
    data_pos = data_pos_all(tv);
    min_pos = min(data_pos)-100;
    max_pos = max(data_pos)+100;
	
	max_rc = max(data_rc);
	rc_step = floor(lambda);
    
    indx1 = find(chr_seg==chr);
    
    %plot
	h=figure(1);
	FontSize = 11;
	set(h,'visible','off');
    clf;
    marker_size = 3;
    %---plot CN---
    subplot(3,1,1);
    hold on
    set(gca,'YGrid','on');
    set(gca,'FontSize',FontSize);
    axis ([min_pos max_pos -0.1 7.1]);
%     set(gca,'XTick',[]);
    set(gca,'YTick',[0:1:7],'Box','on')
    set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
    ylabel('Copy Number','FontSize',FontSize+2);
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        mCN = mcn_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if isempty(indx)
            continue;
        end
        if CN > 7
            CN = 7;
        end
        
        e_indx = length(indx);
        s_indx = 1;
        for k = 2:length(indx)
            if data_pos(indx(k))-data_pos(indx(k-1)) > 10000
                e_indx = k-1;
                s_indx = k;
                break;
            end
        end
        line_style = 'r-';
        plot([data_pos(indx(1)) data_pos(indx(e_indx))],[CN CN],line_style,'LineWidth',2.5);
        plot([data_pos(indx(s_indx)) data_pos(indx(end))],[CN CN],line_style,'LineWidth',2.5);
        line_style = 'b-';
        if CN == mCN
            plot([data_pos(indx(1)) data_pos(indx(e_indx))],[mCN-0.25 mCN-0.25],line_style,'LineWidth',2.5);
            plot([data_pos(indx(s_indx)) data_pos(indx(end))],[mCN-0.25 mCN-0.25],line_style,'LineWidth',2.5);
        else
            plot([data_pos(indx(1)) data_pos(indx(e_indx))],[mCN mCN],line_style,'LineWidth',2.5);
            plot([data_pos(indx(s_indx)) data_pos(indx(end))],[mCN mCN],line_style,'LineWidth',2.5);
        end
    end
    %replace '_' with '-' in barcode to display it correctly
	barcode_m = barcode;
    tmp = strfind(barcode_m,'_');
    barcode_m(tmp) = '-';
    set(gca,'XTick',[])
    title (['Chromosome ' num2str(chr) ', ' barcode_m],'FontSize',FontSize+2)
    
    %---plot BAF---
    subplot(3,1,2);
    set(gca,'FontSize',FontSize);
    hold on
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        if CN == 0
            Muc = 0.5;
        else
            Muc = mcn_seg(j)/cn_seg(j);
        end
        tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
        if sum(tv) == 0
            continue;
        end
        if CN < 2
            k = 1; % Del
        else
            if Muc == 1
                k = 3; % LOH
            else
                k = 2; % Het
            end
        end
        plot(data_pos(tv),data_baf(tv),'.','MarkerSize',marker_size, 'Color', baf_colors(k,:));
    end
%     plot(data_pos,data_baf,'b.', 'MarkerSize',marker_size)
    for j = 0:0.25:1
        plot([data_pos(1) data_pos(end)],[j j],'k-','LineWidth',0.5)
    end
    %plot expected BAF mean values
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        mCN = mcn_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if CN == 0 || mCN/CN == 0.5
            baf_mean = 0.5;
        else
            tv = cns == CN & mcns == mCN;
            if sum(tv) == 0
                baf_mean = mCN/CN;
            else
                baf_mean = mu(tv);
            end
        end
        
        e_indx = length(indx);
        s_indx = 1;
        for k = 2:length(indx)
            if data_pos(indx(k))-data_pos(indx(k-1)) > 10000
                e_indx = k-1;
                s_indx = k;
                break;
            end
        end
        plot([data_pos(indx(1)) data_pos(indx(e_indx))],[baf_mean baf_mean],'k-','LineWidth',1.2);
        plot([data_pos(indx(s_indx)) data_pos(indx(end))],[baf_mean baf_mean],'k-','LineWidth',1.2);
        plot([data_pos(indx(1)) data_pos(indx(e_indx))],[1-baf_mean 1-baf_mean],'k-','LineWidth',1.2);
        plot([data_pos(indx(s_indx)) data_pos(indx(end))],[1-baf_mean 1-baf_mean],'k-','LineWidth',1.2);
    end
    ylabel('B Allele Ratio','FontSize',FontSize+2);
    axis ([min_pos max_pos -0.03 1.03])
%     set(gca,'XTick',[]);
    %set(gca,'YTick',[0 0.5 1],'Box','on')
	set(gca,'YTick',0:0.25:1,'Box','on')
    set(gca,'XTick',[])
    
	%---plot read counts---
    subplot(3,1,3);
    set(gca,'FontSize',FontSize);
    hold on
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
        if sum(tv) == 0
            continue;
        end
        if CN < 1
            k = 1;
        else
            k = CN+1;
        end
        if k > 4
            k = 4;
        end
        plot(data_pos(tv),data_rc(tv),'.','MarkerSize',marker_size, 'Color', rc_colors(k,:));
    end
    % plot expected RC mean values
    for j = reshape(indx1,1,[])
        CN = cn_seg(j);
        indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
        if CN == 0
            CN = 0.01;
        end
        rc_mean = lambda*CN/2;
        
        e_indx = length(indx);
        s_indx = 1;
        for k = 2:length(indx)
            if data_pos(indx(k))-data_pos(indx(k-1)) > 10000
                e_indx = k-1;
                s_indx = k;
                break;
            end
        end
        plot([data_pos(indx(1)) data_pos(indx(e_indx))],[rc_mean rc_mean],'k-','LineWidth',1.2);
        plot([data_pos(indx(s_indx)) data_pos(indx(end))],[rc_mean rc_mean],'k-','LineWidth',1.2);
    end
    
    ylabel('Read Counts','FontSize',FontSize+2);
	set(gca,'YGrid','on','Box','on');
	set(gca,'YTick',0:rc_step:max_rc)
%     axis ([min_pos max_pos min_rd max_rd])
    axis ([min_pos max_pos 0 Inf])
%     set(gca,'XTick',[])
%     axis([-Inf Inf -Inf Inf])

%     subplot(4,1,4);
%     set(gca,'FontSize',FontSize);
%     hold on
%     for j = reshape(indx1,1,[])
%         line_style = 'r-';
%         indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
%         if isempty(indx)
%             continue;
%         end
%         plot([data_pos(indx(1)) data_pos(indx(end))],[score_seg(j) score_seg(j)], line_style, 'LineWidth',1.5);
%     end
%     axis ([min_pos max_pos -10 110])
% %     set(gca,'XTick',[]);
%     set(gca,'YTick',[0:50:100],'Box','on')
%     ylabel('Score','FontSize',FontSize+2);

    %save figure
%     figpath = [plotsdir '\Chr_' num2str(i) '_' barcode];
    figpath = [plotsdir '/Chr_' num2str(chr) '_' barcode '.png'];
    %eval(['print -djpeg -r600 ' figpath ])
    eval(['print -dpng -r400 ' figpath ])
    %print('-dpng','-r400',figpath)
    %saveas(h,figpath,'png');
else
	chromosomes = intersect(unique(data_chr_all),1:24);
	stepsize_ds = 3;
	max_pos = zeros(1,length(chromosomes));
	max_rc = max(data_rc_all);
	min_rc = min(data_rc_all);
	rc_step = floor(lambda);

	for i = 1:length(chromosomes)
		tv1 = data_chr_all == chromosomes(i);
		max_pos(i) = max(data_pos_all(tv1));
	end

	ratio = max_pos/sum(max_pos);
	xtick = cumsum([0 ratio(1:end-1)])+ratio/2;

	FontSize = 17;
	h=figure(1);
	set(h,'visible','off');
	set(gcf,'PaperUnits','inches','PaperPosition',[0 0 13 8])
	clf;
	line_style = 'k-';
	LineWidth = 0.5;
	MarkerSize = 4;
	subplot(3,1,1);
	hold on
	set(gca,'YGrid','on');
	set(gca,'FontSize',FontSize);
	set(gca,'YTick',[0:1:7],'Box','on');
	set(gca,'YTickLabel',{'0','1','2','3','4','5','6','>=7'});
	ylabel('Copy Number');
	pre_x = 0;
	chr_epos = zeros(length(chromosomes),1);
	for i = 1:length(chromosomes)
		tv = data_chr_all == chromosomes(i);
		data_pos = data_pos_all(tv);
		indx = 1:stepsize_ds:length(data_pos);
		data_pos = data_pos(indx);
		x = data_pos*ratio(i)/max_pos(i)+pre_x;
		indx1 = find(chr_seg == chromosomes(i));
		for j = reshape(indx1,1,[])
			CN = cn_seg(j);
			mCN = mcn_seg(j);
			indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
			if isempty(indx)
				continue;
			end
			if CN > 7
				CN = 7;
			end
			plot([x(indx(1)) x(indx(end))],[CN CN], 'r-', 'LineWidth',2.5);
			if CN == mCN
				plot([x(indx(1)) x(indx(end))],[mCN-0.25 mCN-0.25], 'b-', 'LineWidth',2.5);
			else
				plot([x(indx(1)) x(indx(end))],[mCN mCN], 'b-', 'LineWidth',2.5);
			end
		end
		chr_epos(i) = max(x);
		pre_x = pre_x+ratio(i);  
	end
	for i = 1:length(chromosomes)-1
		plot([chr_epos(i) chr_epos(i)], [-0.1 7.1], line_style, 'LineWidth',LineWidth)
	end
	set(gca,'XTick',[])
	axis([0 1 -0.1 7.1]);
	barcode_m = barcode;
	tmp = strfind(barcode_m,'_');
	barcode_m(tmp) = '-';
	title(barcode_m);

	%B allele ratio
	subplot(3,1,2);
	set(gca,'FontSize',FontSize);
	hold on
	pre_x = 0;
	for i = 1:length(chromosomes)
		tv = data_chr_all == chromosomes(i);
		data_pos = data_pos_all(tv);
		data_bd = data_bd_all(tv);
		data_td = data_td_all(tv);
		data_baf = data_bd./data_td;
		indx = 1:stepsize_ds:length(data_pos);
		data_pos = data_pos(indx);
		data_baf = data_baf(indx);
		x = data_pos*ratio(i)/max_pos(i)+pre_x;
		indx1 = find(chr_seg == chromosomes(i));
		for j = reshape(indx1,1,[])
			CN = cn_seg(j);
			if CN == 0
				Muc = 0.5;
			else
				Muc = mcn_seg(j)/cn_seg(j);
			end
			tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
			if sum(tv) == 0
				continue;
			end
			if CN < 2
				k = 1; % Del
			else
				if Muc == 1
					k = 3; % LOH
				else
					k = 2; % Het
				end
			end
			plot(x(tv),data_baf(tv),'.','MarkerSize',MarkerSize, 'Color', baf_colors(k,:));
		end
		for j = 0:0.25:1
			plot([x(1) x(end)],[j j],'k-','LineWidth',0.5)
		end
		%plot expected BAF mean values
		for j = reshape(indx1,1,[])
			CN = cn_seg(j);
			mCN = mcn_seg(j);
			indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
			if isempty(indx)
				continue;
			end
			if CN == 0 || mCN/CN == 0.5
				baf_mean = 0.5;
			else
				tv = cns == CN & mcns == mCN;
				if sum(tv) == 0
					baf_mean = mCN/CN;
				else
					baf_mean = mu(tv);
				end
			end
			plot([x(indx(1)) x(indx(end))],[baf_mean baf_mean],'k-','LineWidth',1.5);
			plot([x(indx(1)) x(indx(end))],[1-baf_mean 1-baf_mean],'k-','LineWidth',1.5);
		end
		pre_x = pre_x+ratio(i);  
	end
	for i = 1:length(chromosomes)-1
		plot([chr_epos(i) chr_epos(i)], [-0.03 1.03], line_style, 'LineWidth',LineWidth)
	end
	ylabel('B Allele Ratio');
	axis([0 1 -0.03 1.03])
	%set(gca,'YTick',[0 0.5 1],'Box','on')
	set(gca,'YTick',0:0.25:1,'Box','on')
	set(gca,'XTick',[])

	subplot(3,1,3);
	set(gca,'FontSize',FontSize);
	hold on
	pre_x = 0;
	for i = 1:length(chromosomes)
		tv = data_chr_all == chromosomes(i);
		data_pos = data_pos_all(tv);
		data_rc = data_rc_all(tv);
		indx = 1:stepsize_ds:length(data_pos);
		data_pos = data_pos(indx);
		data_rc = data_rc(indx);
		x = data_pos*ratio(i)/max_pos(i)+pre_x;
		indx1 = find(chr_seg == chromosomes(i));
		for j = reshape(indx1,1,[])
			CN = cn_seg(j);
			tv = data_pos >= pstart_seg(j) & data_pos <= pend_seg(j);
			if sum(tv) == 0
				continue;
			end
			if CN < 1
				k = 1;
			else
				k = CN+1;
			end
			if k > 4
				k = 4;
			end
			plot(x(tv),data_rc(tv),'.','MarkerSize',MarkerSize, 'Color', rc_colors(k,:));
		end
		% plot expected LCR mean values
		for j = reshape(indx1,1,[])
			CN = cn_seg(j);
			if CN == 0
				CN = 0.01;
			end
			rc_mean = lambda*CN/2;
			indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
			if isempty(indx)
				continue;
			end
			plot([x(indx(1)) x(indx(end))],[rc_mean rc_mean],'k-','LineWidth',1.5);
		end
		pre_x = pre_x+ratio(i);  
	end
	for i = 1:length(chromosomes)-1
		plot([chr_epos(i) chr_epos(i)], [0 max_rc], line_style, 'LineWidth',LineWidth)
	end
	ylabel('Read Counts');
	set(gca,'YGrid','on','Box','on');
	axis([0 1 min_rc max_rc])
	% axis([0 1 -Inf Inf])
	set(gca,'XTick',[])
	set(gca,'YTick',0:rc_step:max_rc);
	xlabel('Chromosome');

	% subplot(4,1,4);
	% set(gca,'FontSize',FontSize);
	% hold on
	% pre_x = 0;
	% for i = 1:length(chromosomes)
	%     tv = data_chr_all == chromosomes(i);
	%     data_pos = data_pos_all(tv);
	%     indx = 1:stepsize_ds:length(data_pos);
	%     data_pos = data_pos(indx);
	%     x = data_pos*ratio(i)/max_pos(i)+pre_x;
	%     indx1 = find(chr_seg == chromosomes(i));
	%     for j = reshape(indx1,1,[])
	%         indx = find(data_pos >= pstart_seg(j) & data_pos <= pend_seg(j));
	%         if isempty(indx)
	%             continue;
	%         end
	%         plot([x(indx(1)) x(indx(end))],[score_seg(j) score_seg(j)], 'r-', 'LineWidth',2.5);
	%     end
	%     pre_x = pre_x+ratio(i);  
	% end
	% for i = 1:length(chromosomes)-1
	%     plot([chr_epos(i) chr_epos(i)], [0 101], line_style, 'LineWidth',LineWidth)
	% end
	% axis ([0 1 0 101])
	% % set(gca,'XTick',[]);
	% set(gca,'YTick',[0:50:100],'Box','on')
	% ylabel('Score');
	% xlabel('Chromosome');

	indx = 1:2:length(chromosomes);
	set(gca,'XTick',xtick(indx));
	set(gca,'XTickLabel',mat2cell(chromosomes(indx)',1,length(indx)));

	%save figure
	figpath = [plotsdir '/' barcode '.png'];
	eval(['print -dpng -r400 ' figpath ])
end

