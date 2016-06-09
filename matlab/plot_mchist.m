%%MATLAB
function out = plot_results(args)

    %close all;

    filename = args.arg1;
    mode = args.arg2;
    
	split=strread(filename,'%s','delimiter','.');
	fname_noext = split{1};

	if strcmp('subplot', mode)
		plotformat = 'subplot';
    elseif strcmp('resolution', mode)
        plotformat = 'resolution';
    else
        plotformat = 'individual';
    end

    %%%%%%%%%% load Monte Carlo results %%%%%%%%%%

    fid1 = fopen(filename);
    sfit = textscan(fid1, '%f64%f64%f64%f64%f64%f64', 'HeaderLines',1);
    fclose(fid1);
    
    
    fid2 = fopen('stimulus.dat');
    parameters = textscan(fid2, '%f64%f64%f64%f64%f64%f64%d%f64%f64%f64');
    fclose(fid2);

    samples = length(sfit{1,1});
%    samples = 100;
    
    estimated.amplitude = sfit{1,1}(1:samples);
    estimated.frequency = sfit{1,2}(1:samples);
    estimated.phase = sfit{1,3}(1:samples);
    estimated.dc = sfit{1,4}(1:samples);
    estimated.erms = sfit{1,5}(1:samples);
    estimated.rtime = sfit{1,6}(1:samples);

    stimulus.amplitude = parameters{1,1}(1:samples);
    stimulus.frequency = parameters{1,2}(1:samples);
    stimulus.phase = parameters{1,3}(1:samples);
    stimulus.dc = parameters{1,4}(1:samples);
    stimulus.fs = parameters{1,5}(1:samples);     % sampling frequency
    stimulus.N = parameters{1,6}(1:samples);      % number of samples
    stimulus.qbits = parameters{1,7}(1:samples);  % number of quantizer bits
    stimulus.Q = parameters{1,8}(1:samples);      % ideal quantizer step size
    stimulus.vmax = parameters{1,9}(1:samples);
    stimulus.vmin = parameters{1,10}(1:samples);

    %%%%%%%%%% prepare histogram data %%%%%%%%%%
    Delta.A=(stimulus.amplitude-estimated.amplitude.*stimulus.Q)./stimulus.Q;
    temp = Delta.A;
    Delta_mean.A = mean(temp);
    Delta_std.A = std(temp);

    Delta.f=(stimulus.frequency-estimated.frequency)./stimulus.frequency;
    temp = Delta.f;
    Delta_mean.f = mean(temp);
    Delta_std.f = std(temp);

    Delta.phi=(stimulus.phase-estimated.phase)./stimulus.phase;
    temp = Delta.phi;
    Delta_mean.phi = mean(temp);
    Delta_std.phi = std(temp);

    Delta.dc=(stimulus.dc-(estimated.dc.*stimulus.Q+stimulus.vmin))./stimulus.Q;
    temp = Delta.dc;
    Delta_mean.dc = mean(temp);
    Delta_std.dc = std(temp);

    Delta.rtime=(estimated.rtime)./stimulus.N;
    temp = Delta.rtime;
    Delta_mean.rtime = mean(temp);
    Delta_std.rtime = std(temp);

    %%%%%%%%%% plot histograms %%%%%%%%%%

    if strcmp('individual',plotformat)

        % plot estimation error of the amplitude
        h0 = figure;
        hist(Delta.A,length(Delta.A)/3);
        hp0 = findobj(gca,'Type','patch');
        set(hp0,'FaceColor','black','EdgeColor','none');
        grid on;
        hc0 = get(h0,'children');
        set(hc0,'fontsize',16);
        legend(sprintf('mean: %g\nstd: %g',Delta_mean.A, Delta_std.A),'Location','NE');
        title('Estimation Error of the Amplitude');
        xlabel('(A-A\_est)/Q');
        ylabel('Frequency');
        qfigset;
        print('-deps', sprintf('%s_amp_hist.eps', fname_noext));
        print('-dpng', sprintf('%s_amp_hist.png', fname_noext));


        % plot estimation error of the frequency
        h1 = figure;
        hist(Delta.f,length(Delta.f)/3);
        hp1 = findobj(gca,'Type','patch');
        set(hp1,'FaceColor','black','EdgeColor','none');
        grid on;
        hc1 = get(h1,'children');
        set(hc1,'fontsize',16);
        legend(sprintf('mean: %g\nstd: %g',Delta_mean.f, Delta_std.f),'Location','NE');
        title('Estimation Error of the Frequency');
        xlabel('(f-f\_est)/f');
        ylabel('Frequency');
        qfigset;
        print('-deps', sprintf('%s_freq_hist.eps', fname_noext));
        print('-dpng', sprintf('%s_freq_hist.png', fname_noext));

 
        % plot estimation error of the phase
        h2 = figure;
        hist(Delta.phi,length(Delta.phi)/3);
        hp2 = findobj(gca,'Type','patch');
        set(hp2,'FaceColor','black','EdgeColor','none');
        grid on;
        hc2 = get(h2,'children');
        set(hc2,'fontsize',16);
        legend(sprintf('mean: %g\nstd: %g',Delta_mean.phi, Delta_std.phi),'Location','NE');
        title('Estimation Error of the Phase');
        xlabel('(phase-phase\_est)/phase');
        ylabel('Frequency');
        qfigset;
        print('-deps', sprintf('%s_phase_hist.eps', fname_noext));
        print('-dpng', sprintf('%s_phase_hist.png', fname_noext));


        % plot estimation error of the dc
        h4 = figure;
        hist(Delta.dc,length(Delta.dc)/3);
        hp4 = findobj(gca,'Type','patch');
        set(hp4,'FaceColor','black','EdgeColor','none');
        grid on;
        hc4 = get(h4,'children');
        set(hc4,'fontsize',16);
        legend(sprintf('mean: %g\nstd: %g',Delta_mean.dc, Delta_std.dc),'Location','NE');
        title('Estimation Error of the dc');
        xlabel('(dc-dc\_est)/Q');
        ylabel('Frequency');
        qfigset;
        print('-deps', sprintf('%s_dc_hist.eps', fname_noext));
        print('-dpng', sprintf('%s_dc_hist.png', fname_noext));
 

        % plot runtime histogram
        h4 = figure;
        hist(Delta.rtime,length(Delta.rtime)/3);
        hp4 = findobj(gca,'Type','patch');
        set(hp4,'FaceColor','black','EdgeColor','none');
        grid on;
        hc4 = get(h4,'children');
        set(hc4,'fontsize',16);
        legend(sprintf('mean: %g\nstd: %g',Delta_mean.rtime, Delta_std.rtime),'Location','NE');
        title('Normalized Runtime');
        xlabel('(rtime)/N');
        ylabel('Frequency');
        qfigset;
        print('-deps', sprintf('%s_rtime_hist.eps', fname_noext));
        print('-dpng', sprintf('%s_rtime_hist.png', fname_noext));

 
    elseif strcmp('subplot',plotformat)

        fig = figure;
        % plot estimation error of the amplitude
        [n1,x1]=hist(Delta.A,length(Delta.A)/3);
        subplot(2,2,1), bar(x1,n1);
        hp5 = findobj(gca,'Type','patch');
        set(hp5,'FaceColor','black','EdgeColor','none');
        grid on;
        legend(hp5,sprintf('mean: %g\nstd: %g',Delta_mean.A, Delta_std.A),'Location','NE');
        title('Estimation Error of the Amplitude');
        xlabel('(A-A\_est)/Q');
        ylabel('Frequency');

        % plot estimation error of the frequency
        [n2,x2]=hist(Delta.f,length(Delta.f)/3);
        subplot(2,2,2), bar(x2,n2);
        hp6 = findobj(gca,'Type','patch');
        set(hp6,'FaceColor','black','EdgeColor','none');
        grid on;
        legend(hp6, sprintf('mean: %g\nstd: %g',Delta_mean.f, Delta_std.f),'Location','NE');
        title('Estimation Error of the Frequency');
        xlabel('(f-f\_est)/f');
        ylabel('Frequency');

        % plot estimation error of the phase
        [n3,x3]= hist(Delta.phi,length(Delta.phi)/3);
        subplot(2,2,3), bar(x3,n3);
        hp7 = findobj(gca,'Type','patch');
        set(hp7,'FaceColor','black','EdgeColor','none');
        grid on;
        legend(hp7, sprintf('mean: %g\nstd: %g',Delta_mean.phi, Delta_std.phi),'Location','NE');
        title('Estimation Error of the Phase');
        xlabel('(phase-phase\_est)/phase');
        ylabel('Frequency');

        % plot estimation error of the dc
        [n4,x4]=hist(Delta.dc,length(Delta.dc)/3);
        subplot(2,2,4), bar(x4,n4);
        hp8 = findobj(gca,'Type','patch');
        set(hp8,'FaceColor','black','EdgeColor','none');
        grid on;
        legend(hp8, sprintf('mean: %g\nstd: %g',Delta_mean.dc, Delta_std.dc),'Location','NE');
        title('Estimation Error of the dc');
        xlabel('(dc-dc\_est)/Q');
        ylabel('Frequency');

        % resize figure using the golden ratio
        gold=(sqrt(5)-1)/2; %0.6180
        pos=get(fig, 'OuterPosition');
        pos(1)=200;
        pos(2)=150;
        pos(4)=850;
        pos(3)=round(pos(4)/gold);
        set(fig, 'Units', 'pixels', 'OuterPosition', pos);
        % to avoid resizing the figure during print
        set(fig, 'PaperPositionMode', 'auto');

        print('-deps', sprintf('%s_hist.eps', fname_noext));
        print('-dpng', sprintf('%s_hist.png', fname_noext));

        
        fig = figure;
        % plot estimation error of the amplitude
        subplot(2,2,1), plot(Delta.A,'.-');
        grid on;
        title(sprintf('Estimation Error of the Amplitude'));
        xlabel('Monte Carlo simulations');
        ylabel('(A-A\_est)/Q');
        
        % plot estimation error of the frequency
        subplot(2,2,2), plot(Delta.f,'.-');
        grid on;
        title(sprintf('Estimation Error of the Frequency'));
        xlabel('Monte Carlo simulations');
        ylabel('(f-f\_est)/f');        
        
        % plot estimation error of the phase
        subplot(2,2,3), plot(Delta.phi,'.-');
        grid on;
        title(sprintf('Estimation Error of the Phase'));
        xlabel('Monte Carlo simulations');
        ylabel('(phase-phase\_est)/phase');        
        
        % plot estimation error of the dc
        subplot(2,2,4), plot(Delta.dc,'.-');
        grid on;
        title(sprintf('Estimation Error of the DC'));
        xlabel('Monte Carlo simulations');
        ylabel('(dc-dc\_est)/Q');
        
        % resize figure using the golden ratio
        gold=(sqrt(5)-1)/2; %0.6180
        pos=get(fig, 'OuterPosition');
        pos(1)=200;
        pos(2)=150;
        pos(4)=850;
        pos(3)=round(pos(4)/gold);
        set(fig, 'Units', 'pixels', 'OuterPosition', pos);
        % to avoid resizing the figure during print
        set(fig, 'PaperPositionMode', 'auto');
        
        print('-deps', sprintf('%s_seq.eps', fname_noext));
        print('-dpng', sprintf('%s_seq.png', fname_noext));
        
        
    elseif strcmp('resolution',plotformat)

        fig = figure;
        % plot estimation error of the amplitude
        subplot(2,2,1), plot(stimulus.vmax, Delta.A,'.-');
        grid on;
        title('Estimation Error of the Amplitude vs Vmax');
        xlabel('Vmax');
        ylabel('(A-A\_est)/Q');

        % plot estimation error of the frequency
        subplot(2,2,2), plot(stimulus.vmax, Delta.f,'.-');
        grid on;
        title('Estimation Error of the Frequency vs Vmax');
        xlabel('Vmax');
        ylabel('(f-f\_est)/f');

        % plot estimation error of the phase
        subplot(2,2,3), plot(stimulus.vmax, Delta.phi,'.-');
        grid on;
        title('Estimation Error of the Phase vs Vmax');
        xlabel('Vmax');
        ylabel('(phase-phase\_est)/phase');

        % plot estimation error of the dc
        subplot(2,2,4), plot(stimulus.vmax, Delta.dc,'.-');
        grid on;
        title('Estimation Error of the dc vs Vmax');
        xlabel('Vmax');
        ylabel('(dc-dc\_est)/Q');

        % resize figure using the golden ratio
        gold=(sqrt(5)-1)/2; %0.6180
        pos=get(fig, 'OuterPosition');
        pos(1)=200;
        pos(2)=150;
        pos(4)=850;
        pos(3)=round(pos(4)/gold);
        set(fig, 'Units', 'pixels', 'OuterPosition', pos);
        % to avoid resizing the figure during print
        set(fig, 'PaperPositionMode', 'auto');

        print('-deps', sprintf('%s_resvsVmax.eps', fname_noext));
        print('-dpng', sprintf('%s_resvsVmax.png', fname_noext));

    end

    disp('done');
end
