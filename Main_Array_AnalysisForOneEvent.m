clc
clear
close all
%%% the sac headers must include the true back-azimuth, check the baz
%%% header
addpath(genpath('./src/'))
! rm -rf ./Figures/*
file_path = './Waveform/';
wavdir = 'ev.list';
% list   = 'list';
%% example3
EVID   ='wenchuan';
% afterB = 450; window =20 ;           %(s)
afterB = 454; window =7;
% fmin   = 0.01; fmax = 1; fstep =0.01;  %(Hz)
fmin   = 0.1; fmax = 0.6; fstep =0.01;  %(Hz)
%% example2
% EVID='2016-12-08-mwb65-off-coast-of-northern-california';
% afterB = 55;
% window =20;
% fmin=5;
% fmax=15;
% fstep=0.01;
%% example3
% EVID='2016-12-14-mw50-northern-california';
% EVID='ID_rename_ev_12_211';
% afterB = 3.1;
% window =0.58;
% fmin=5;
% fmax=15;
% fstep=0.01;
slowmin = 0.0;
slowmax = 0.2;

bazmax = 360; bazmin = 0;
sgrids    = 51;
baz_grids = 361;
dslow     = (slowmax-slowmin)/(sgrids-1);
% dslow     = slowmax/(sgrids-1);
dbaz   = (bazmax-bazmin)/(baz_grids-1);
slow   = (0:1:sgrids-1);
baz    = (0:1:baz_grids-1);
% slow   = slow*dslow;
slow   = slowmin + slow*dslow;
baz    = bazmin+ baz*dbaz;
tgrids = baz_grids*sgrids;
% return
%% get the record
[data,STA,wavedata] = process_data(wavdir,EVID,afterB,window);
%%
% figure
% spectrogram(wavedata(9,1).dataE(:),128,100,512,100,'yaxis')
% ylim([0 1])
% xlim([0 1])
% return
sta_xy(:,1) = STA.x/1000; sta_xy(:,2) = STA.y/1000;  % m convert to km
swd         = data;                                  % dimension is (samples,nsta)
staName     = STA.name;
delta       = STA.delta(1);
nsta        = size(sta_xy,1);
%% construct spectra  measurements
[Amp,freqs,nfft] = fft_matrix(swd,delta);
% figure; % plot spectrum from function fft_matrix
% set(gcf,'PaperPositionMode','auto')
% set(gcf,'units','centimeters')
% set(gcf,'pos',[8 12 20 15])
% semilogy(freqs,abs(Amp(:,STA.refindex)));
% loglog(freqs,abs(Amp(:,STA.refindex)),'-k','linewidth',2);
% hold on
% % vert=[0.3 10e-3;0.3 10e7;0.6 10e7;0.6 10e-3];
% % face =[1,2,3,4];
% % patch('Faces',face,'Vertices',vert,...
% %     'EdgeColor','none','FaceColor','k','LineWidth',0.1);
% xlabel('Frequency (Hz)')
% ylabel('Ampitude')
% ylim([10e-3 10e7])
% % set(gca,'ytick',[10e-3 10e0  10e4  10e7])
% set(gca,'linewidth',2,'fontsize',22)
% alpha(0.25)
% print(gcf,'-depsc2','-r600','FrequencyBand.eps')
% return
%% modeling time table with slowness and backazimuth range
[delay_table]=get_delaytable(slow,baz,sta_xy(:,1),sta_xy(:,2)) ;
%% select frequency steps according to the starting fmin and fmax
[selfreq,freq_index,fsstep,nfreq]=get_bandwidth(freqs,fmin,fmax,fstep);

alpha   = 1.1 ;%0.05 0.1  1.0 1.1 1.3];
% alpha = sqrt(nsta)*0.0001;
retol   = 0.01;
x       = zeros(2*tgrids,1) ;
Xw      = zeros(2*tgrids,nfreq) ;
tgrids2 = 2*tgrids;
%% loop all frequency steps
for i = 1:nfreq
    fc = selfreq(i);
    fprintf('Inverse frequency  %5.3f Hz\n',fc)
    
    Bw  = [real(Amp(freq_index(1)+i,:))';imag(Amp(freq_index(1)+i,:))'];
    [G] = get_cs_matrix(nsta,delay_table,tgrids,fc);
    %     % mutiple a random matrix
    %             FAI=randn(4*nsta,2*nsta);
    %             G=FAI*G;
    %             Bw=FAI*Bw;
    
    
    %    OMP Method (is very fast)
    
    x = cs_omp(Bw,G,tgrids2);
    
    %    l1_ls Method (worse)
    %
    %           [x] = l1_ls(G,Bw,alpha,retol);
    
    %     CVX  Method (is very slowly)
    %     if test it, CVX package should be installed in Matlab. 
    
    %     cvx_begin
    %     variable x(2*tgrids)
    %     minimize( norm(G*x-Bw,2)+alpha*norm(x,1) );
    %     cvx_end
    %
    %     BCS   Method  (may be MORE robust than OMP)
    
    %     initsigma2 = std(Bw)^2/1e2;
    %     [weights,used,sigma2,errbars] = BCS_fast_rvm(G,Bw,initsigma2,1e-8);
    %     x = zeros(tgrids2,1);
    %     x(used) = weights;
    
    Xw(:,i) = x(:);
end
%% store all results
Estored = zeros(tgrids,nfreq);
for i = 1:nfreq
    re_X = Xw(1:tgrids,i);
    im_X = Xw(tgrids+1:end,i);
    AMP = complex(re_X,im_X);
    Estored(:,i) = sqrt(re_X.^2+im_X.^2);
end
%% all results are stored in E(tgrids,nfreqs) array using variant frequency
%  we choose to show each frequency result or show result from the
%  reduction of all energy of each frequency

%%
show_results = 'mutiple';
% show_results = 'reduction';
plot_flag ='true';
STA.baz(1)=272.2;   % should not be in there if the sac headers include true baz header.
baz_true=STA.baz(1);
if (strcmp(show_results,'mutiple'))  % show each frequency result
    for i=1:nfreq
        fc = selfreq(i);
        fprintf('\nResult from frequency  %5.3f Hz\n',fc)
[fk_energy,baz_fk,slo_fk]=fk_analysis_polar_singlefrequency(Amp,freqs,fc,baz,slow,sta_xy);        
        E = Estored(:,i);
%         figure(56)
%         clf
%         plot(E)
%         pause(1)
        showresults(baz,slow,E,sgrids,baz_grids,fk_energy,baz_true,EVID,fc,plot_flag);
    end
elseif(strcmp(show_results,'reduction')) % show reduction result
    E = zeros(tgrids,1);
    for i=1:nfreq
        E = E + Estored(:,i);
    end
    % FK result
    [fk_energy,~,~]=fk_analysis_polar(Amp,freqs,fmin,fmax,baz,slow,sta_xy);
    [baz_cs,slo_cs,baz_fk,slo_fk]= showresults(baz,slow,E,sgrids,baz_grids,fk_energy,baz_true,EVID,1,plot_flag);
end
