%% MATLAB
function out = plot_deltahist(args)


% This script calculates the differences of the estimation
% results of two algorithms from the Monte Carlo simulation,
% also plots the delta sequence and the histogram.
%
%       out = plot_deltahist(args)
%        
%       Input arguments:
%         args.arg1: the name of the first file to process
%         args.arg2: the name of the second file to process
%
%       This program is public domain. It is available through
%       https://github.com/asarhegyi/

% $Id:  $
% Copyright (c) 2015-2016 by Attila Sarhegyi
% All rights reserved.


    %close all;

    filename1 = args.arg1;
    filename2 = args.arg2;
    
    split1=strread(filename1,'%s','delimiter','.');
    fname1_noext = split1{1};
    title1a = strrep(fname1_noext, '_', ' ');
    title1b = strrep(title1a, 'cpp', 'C++');
    title1c = strrep(title1b, 'matlab', 'MATLAB');
    
    split=strread(filename2,'%s','delimiter','.');
    fname2_noext = split{1};
    title2a = strrep(fname2_noext, '_', ' ');
    title2b = strrep(title2a, 'cpp', 'C++');
    title2c = strrep(title2b, 'matlab', 'MATLAB');

    %%%%%%%%%% load Monte Carlo results %%%%%%%%%%

    fida = fopen(filename1);
    sfit_a = textscan(fida, '%f64%f64%f64%f64%f64%f64', 'HeaderLines',1);
    fclose(fida);

    fidb = fopen(filename2);
    sfit_b = textscan(fidb, '%f64%f64%f64%f64%f64%f64', 'HeaderLines',1);
    fclose(fidb);

    est_a.amplitude = sfit_a{1,1};
    est_a.frequency = sfit_a{1,2};
    est_a.phase = sfit_a{1,3};
    est_a.dc = sfit_a{1,4};
%    est_a.erms = sfit_a{1,5};
%    est_a.rtime = sfit_a{1,6};

    est_b.amplitude = sfit_b{1,1};
    est_b.frequency = sfit_b{1,2};
    est_b.phase = sfit_b{1,3};
    est_b.dc = sfit_b{1,4};
%    est_b.erms = sfit_b{1,5};
%    est_b.rtime = sfit_b{1,6};

    %%%%%%%%%% prepare histogram data %%%%%%%%%%

    delta_cpp.A=(est_a.amplitude - est_b.amplitude)./est_a.amplitude;
    temp = delta_cpp.A;
    delta_cpp_mean.A = mean(temp);
    delta_cpp_std.A = std(temp);

    delta_cpp.f=(est_a.frequency - est_b.frequency)./est_a.frequency;
    temp = delta_cpp.f;
    delta_cpp_mean.f = mean(temp);
    delta_cpp_std.f = std(temp);

    delta_cpp.phi=(est_a.phase - est_b.phase)./est_a.phase;
    temp = delta_cpp.phi;
    delta_cpp_mean.phi = mean(temp);
    delta_cpp_std.phi = std(temp);

    delta_cpp.dc=(est_a.dc - est_b.dc)./est_a.dc;
    temp = delta_cpp.dc;
    delta_cpp_mean.dc = mean(temp);
    delta_cpp_std.dc = std(temp);

    %%%%%%%%%% plot histograms %%%%%%%%%%

    fig = figure;
    % plot estimation error of the amplitude
    [n1,x1]=hist(delta_cpp.A,length(delta_cpp.A)/3);
    subplot(2,2,1), bar(x1,n1);
    hp5 = findobj(gca,'Type','patch');
    set(hp5,'FaceColor','black','EdgeColor','none');
    grid on;
    legend(hp5,sprintf('mean: %g\nstd: %g',delta_cpp_mean.A, delta_cpp_std.A),'Location','NE');
    title(sprintf('Estimation Error of the Amplitude: %s vs. %s', title1c, title2c));
    xlabel('(A\_mat-A\_cpp)/A\_mat');
    ylabel('Frequency');

    % plot estimation error of the frequency
    [n2,x2]=hist(delta_cpp.f,length(delta_cpp.f)/5);
    subplot(2,2,2), bar(x2,n2);
    hp6 = findobj(gca,'Type','patch');
    set(hp6,'FaceColor','black','EdgeColor','none');
    grid on;
    legend(hp6, sprintf('mean: %g\nstd: %g',delta_cpp_mean.f, delta_cpp_std.f),'Location','NE');
    title(sprintf('Estimation Error of the Frequency: %s vs. %s', title1c, title2c));
    xlabel('(f\_mat-f\_cpp)/f\_mat');
    ylabel('Frequency');

    % plot estimation error of the phase
    [n3,x3]= hist(delta_cpp.phi,length(delta_cpp.phi)/5);
    subplot(2,2,3), bar(x3,n3);
    hp7 = findobj(gca,'Type','patch');
    set(hp7,'FaceColor','black','EdgeColor','none');
    grid on;
    legend(hp7, sprintf('mean: %g\nstd: %g',delta_cpp_mean.phi, delta_cpp_std.phi),'Location','NE');
    title(sprintf('Estimation Error of the Phase: %s vs. %s', title1c, title2c));
    xlabel('(phase\_mat-phase\_cpp)/phase\_mat');
    ylabel('Frequency');

    % plot estimation error of the dc
    [n4,x4]=hist(delta_cpp.dc,length(delta_cpp.dc)/5);
    subplot(2,2,4), bar(x4,n4);
    hp8 = findobj(gca,'Type','patch');
    set(hp8,'FaceColor','black','EdgeColor','none');
    grid on;
    legend(hp8, sprintf('mean: %g\nstd: %g',delta_cpp_mean.dc, delta_cpp_std.dc),'Location','NE');
    title(sprintf('Estimation Error of the DC: %s vs. %s', title1c, title2c));
    xlabel('(dc\_mat-dc\_cpp)/dc\_mat');
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

    print('-deps', sprintf('%s_vs_%s_hist.eps', fname1_noext, fname2_noext));
    print('-dpng', sprintf('%s_vs_%s_hist.png', fname1_noext, fname2_noext));

    %%%%%%%%%% plot sequence %%%%%%%%%%

    fig = figure;
    % plot estimation error of the amplitude
    subplot(2,2,1), plot(delta_cpp.A,'.-');
    grid on;
    title(sprintf('Estimation Error of the Amplitude: %s vs. %s', title1c, title2c));
    xlabel('Monte Carlo simulations');
    ylabel('(A\_mat-A\_cpp)/A\_mat');

    % plot estimation error of the frequency
    subplot(2,2,2), plot(delta_cpp.f,'.-');
    grid on;
    title(sprintf('Estimation Error of the Frequency: %s vs. %s', title1c, title2c));
    xlabel('Monte Carlo simulations');
    ylabel('(f\_mat-f\_cpp)/f\_mat');

    % plot estimation error of the phase
    subplot(2,2,3), plot(delta_cpp.phi,'.-');
    grid on;
    title(sprintf('Estimation Error of the Phase: %s vs. %s', title1c, title2c));
    xlabel('Monte Carlo simulations');
    ylabel('(phase\_mat-phase\_cpp)/phase\_mat');

    % plot estimation error of the dc
    subplot(2,2,4), plot(delta_cpp.dc,'.-');
    grid on;
    title(sprintf('Estimation Error of the DC: %s vs. %s', title1c, title2c));
    xlabel('Monte Carlo simulations');
    ylabel('(dc\_mat-dc\_cpp)/dc\_mat');

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

    print('-deps', sprintf('%s_vs_%s_seq.eps', fname1_noext, fname2_noext));
    print('-dpng', sprintf('%s_vs_%s_seq.png', fname1_noext, fname2_noext));

    disp('done');
end
