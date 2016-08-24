% This script establishes an environment to run the Maximum Likelihood Estimation script
% under development.  

close all;
clear all;
clc;

diary ML_est.log
diary on;

Vmax =  2;
Vmin = -2;

rm.verbose = 2; %run modifier
%rm.verbose = 3; %run modifier
rm.eps_limit = 0.05;
rm.Tl_alg = 'groupped';
%rm.tl_alg = 'full';


MC_run = 500; %Number of Monte Carlo runs
N=8;          %Number of bits

%Initial values for the input sine wave and the noise
freq = 1;
dc_level = 0;
amplitude = 2.021;
Ts = 1/10e3;
N_samples = 8*8000;
sigma = .018;


    rm.noise_model = 'Gauss';
    MC_results = zeros(MC_run,5);

   
    
    for search_model=1:8
    
        switch search_model
            case {1}
                rm.search_alg = 'fminsearch'; rm.sine = 'ABC';    rm.noise = 0;
            case {2}
                rm.search_alg = 'fminsearch'; rm.sine = 'ABC';    rm.noise = 1;
            case {3}
                rm.search_alg = 'fminsearch'; rm.sine = 'AphiDC'; rm.noise = 0;
            case {4}
                rm.search_alg = 'fminsearch'; rm.sine = 'AphiDC'; rm.noise = 1;
            case {5}
                rm.search_alg = 'fminunc';    rm.sine = 'ABC';    rm.noise = 0;
            case {6}
                rm.search_alg = 'fminunc';    rm.sine = 'ABC';    rm.noise = 1;
            case {7}
                rm.search_alg = 'fminunc';    rm.sine = 'AphiDC'; rm.noise = 0;
            case {8}
                rm.search_alg = 'fminunc';    rm.sine = 'AphiDC'; rm.noise = 1;
        end

        log_message = sprintf('Monte Carlo Sim for %s %s %s N%d',rm.noise_model,rm.sine,rm.search_alg,rm.noise);
        disp(log_message);

        for k=1:MC_run

            disp(sprintf('Monte Carlo Iteration #: %d ',k));

            meas_data = mdata;
            meas_data = set_V_max(meas_data,Vmax);
            meas_data = set_V_min(meas_data,Vmin);
            meas_data = set_Nbit(meas_data,N);
            meas_data = set_Ts(meas_data,Ts);
            meas_data = set_sine_freq(meas_data,freq);

            Tl_ideal = get_Tl_of_an_ideal_quantizer(meas_data);

            tic;            %start timer
            toc_temp = toc;
            [bb,cc] = gen_samples(Tl_ideal,sigma,amplitude,freq,dc_level,Ts,N_samples);

            NN = length(bb);
            tt_bb = 0:NN-1;
            tt_bb = tt_bb.';
            bb2 = [bb(1:end/4); bb(3*end/4:end)];
            tt_bb2 = [tt_bb(1:end/4); tt_bb(3*end/4:end)];
            meas_data = set_meas_data(meas_data,bb,tt_bb);

            %%%%%%%%%%%%%%%%%%%%%%%% Plot/Debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if rm.verbose > 2
                figure; plot(tt_bb,bb); grid; hold; plot(20*cc,'r');
            end%if
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [A_est,B_est,C_est,Tl_est,sigma_est] = ML_est_ABCN(meas_data,rm);
            toc_save = toc-toc_temp;

            MC_results(k,:) = [A_est,B_est,C_est,sigma_est,toc_save];

        end
        
        %Save results into a file
        FILEO = sprintf('MC_results_%s_%s_%s_N%d.txt',rm.noise_model,rm.sine,rm.search_alg,rm.noise);
        disp(FILEO);
        save(FILEO,'MC_results','-ASCII'); 
        
        %Calculate the estimation error for the sinewave parameteres
        Sine  = MC_results(:,1:3);
        Noise = MC_results(:,4);
        Time  = MC_results(:,5);
        
        % Input sine wave (noiseless)
        tt = 0:N_samples-1;
        tt = tt*Ts;
        zz = amplitude*sin(2*pi*freq*tt+pi*1.2)+dc_level;

        
        if strcmp(rm.sine,'ABC')
            A = sqrt( Sine(:,1).^2 + Sine(:,2).^2 );
            for k = 1:length(Sine) 
                if      Sine(k,2)>0
                    phi(k)=atan(-Sine(k,1)/Sine(k,2));
                elseif Sine(k,2)<0
                    phi(k)=atan(-Sine(k,1)/Sine(k,2))+pi;
                else %Sine(1)==0
                    if Sine(k,1)<0, phi(k)=pi/2;
                    else phi(k)=3*pi/2;
                    end
                end
            end
            DC = Sine(:,3);
        else
            A   = Sine(:,1);
            phi = Sine(:,2);
            DC  = Sine(:,3);
        end
        
        Delta.A   = (amplitude - A)/amplitude;
        Delta.phi = (pi*1.2 - phi)/(pi*1.2);
        Delta.DC  = (dc_level - DC);
       
        Delta_mean.A   = mean(Delta.A);
        Delta_mean.phi = mean(Delta.phi);
        Delta_mean.DC  = mean(Delta.DC);
        
        Delta_std.A   = std(Delta.A);
        Delta_std.phi = std(Delta.phi);
        Delta_std.DC  = std(Delta.DC);

        close all;
        
        % plot estimation error for Amplitude
        h1 = figure;
        hist(Delta.A,length(Delta.A)/5);
        hp1 = findobj(gca,'Type','patch');
        set(hp1,'FaceColor','black','EdgeColor','none');
        grid on;
        hc1 = get(h1,'children');
        set(hc1,'fontsize',16);
        Title(sprintf('Estimation error of the amplitude\n m=%.2e std=%.2e',Delta_mean.A,Delta_std.A));
        xlabel('(A-A_{est})/A');
        ylabel('Frequency');
        qfigset;
        FILEO = sprintf('MC_results_%s_%s_%s_N%d_Amplitude.eps',rm.noise_model,rm.sine,rm.search_alg,rm.noise);
        disp(FILEO);
        print(FILEO,'-deps');
        
        % plot estimation error for Phase 
        h2 = figure;
        hist(Delta.phi,length(Delta.phi)/5);
        hp2 = findobj(gca,'Type','patch');
        set(hp2,'FaceColor','black','EdgeColor','none');
        grid on;
        hc2 = get(h2,'children');
        set(hc2,'fontsize',16);
        Title(sprintf('Estimation error of the Phase\n m=%.2e std=%.2e',Delta_mean.phi,Delta_std.phi));
        xlabel('(phi-phi_{est})/phi');
        ylabel('Frequency');
        qfigset;
        FILEO = sprintf('MC_results_%s_%s_%s_N%d_Phase.eps',rm.noise_model,rm.sine,rm.search_alg,rm.noise);
        disp(FILEO);
        %print -depsc2 MC_Amp_GvsL.eps;
        print(FILEO,'-deps');
        
        % plot estimation error for DC 
        h3 = figure;
        hist(Delta.DC,length(Delta.DC)/5);
        hp3 = findobj(gca,'Type','patch');
        set(hp3,'FaceColor','black','EdgeColor','none');
        grid on;
        hc3 = get(h3,'children');
        set(hc3,'fontsize',16);
        Title(sprintf('Estimation error of the DC\n m=%.2e std=%.2e',Delta_mean.DC,Delta_std.DC));
        xlabel('(DC-DC_{est})/DC');
        ylabel('Frequency');
        qfigset;
        FILEO = sprintf('MC_results_%s_%s_%s_N%d_DC.eps',rm.noise_model,rm.sine,rm.search_alg,rm.noise);
        disp(FILEO);
        %print -depsc2 MC_Amp_GvsL.eps;
        print(FILEO,'-deps');
       
    end
    
close all;
diary off;