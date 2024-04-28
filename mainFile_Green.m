clc;
clear all;
close all;

%% Parameters initialization
Ind_plotDebug = 0;
Ind_plotting = 1;
fontSize = 10;
lineWidth = 2;
numSamplesPerSecond = 15;
NFFT = 256;
maxFundamentalFreq = 120;   % bpm
minFundamentalFreq = 40;    % bpm
stdDev_HR = 0.05;            % Hz

% fileName = 'data_Z.rtf';
fileName = 'dataCapture_K5_15sps.csv';
% fileName = 'dataCapture_K2.csv';

if strcmp( fileName(end-2:end) , 'csv' )
    %% CSV File read
    rawData = readmatrix(fileName , 'FileType' , 'text' , 'Delimiter' , {','});     % for csv files
    
elseif strcmp( fileName(end-2:end) , 'rtf' )
    %% RTF File read
    rawData = readmatrix(fileName , 'FileType' , 'text' , 'Delimiter' , {',','\'});     % for rtf files
    rawData = rawData( 2:end-1 , 1:3 );

else
    error("Enter a valid file type!")
end

[row , col] = size(rawData);
if col == 4
    vecTime_recorded = rawData(: , 1);
    rawData = rawData(: , 2:4);
end

%% Plots
vecTime = [0: 1/numSamplesPerSecond : (row-1)/numSamplesPerSecond];

figure
plot( vecTime , rawData(:,1) , 'Linewidth', lineWidth )
xlabel( 'Time (s)' );
ylabel( 'Magnitude' );
set( findall( gcf , '-property' , 'FontSize' ) , 'FontSize' , fontSize );
set( gcf , 'color' , [1 1 1] );
legend( "red" , 'FontSize' , fontSize , 'Location' , 'northeast' , 'Box' , 'off' );
grid on;

figure
plot( vecTime , rawData(:,2) , 'Linewidth', lineWidth )
xlabel( 'Time (s)' );
ylabel( 'Magnitude' );
set( findall( gcf , '-property' , 'FontSize' ) , 'FontSize' , fontSize );
set( gcf , 'color' , [1 1 1] );
legend( "Infrared" , 'FontSize' , fontSize , 'Location' , 'northeast' , 'Box' , 'off' );
grid on;

figure
plot( vecTime , rawData(:,3) , 'Linewidth', lineWidth )
xlabel( 'Time (s)' );
ylabel( 'Magnitude' );
set( findall( gcf , '-property' , 'FontSize' ) , 'FontSize' , fontSize );
set( gcf , 'color' , [1 1 1] );
legend( "green" , 'FontSize' , fontSize , 'Location' , 'northeast' , 'Box' , 'off' );
grid on;

%% Algorithm for movement artifacts
% parameters
ts = 1/numSamplesPerSecond;
windowSize = 4;     % s
windowOverlap = 0;  % s
totalTimeLength = floor( (row-1) * ts );

numSamplesPerWindow = ceil( windowSize * numSamplesPerSecond );
windowSize = numSamplesPerWindow * ts;
numSamplesPerOverlap = ceil( windowOverlap * numSamplesPerSecond );
windowOverlap = numSamplesPerOverlap * ts;

numTotalWindows = floor( (totalTimeLength - windowOverlap) / (windowSize - windowOverlap) );
numSamplesNewPerWindow = numSamplesPerWindow - numSamplesPerOverlap;

% elliptic filter
fs = numSamplesPerSecond;
n = 4; 
Rp = 3;
Rs = 50;
Wn = [0.4 4]/(fs/2); 
ftype = 'bandpass'; 
[b,a] = ellip(n,Rp,Rs,Wn,ftype);
if Ind_plotDebug
    freqz(b,a,[],fs);
end

% figures
if Ind_plotDebug
    figure(10)
    hold on;
    title('IR')
    
    figure(11)
    hold on;
    title('Green')
    
    figure(12)
    hold on;
    title('Subtraction')
    
    figure(13)
    hold on;
    title('Green - diff 0')
    
    figure(14)
    hold on;
    title('Green - diff 1')

    figure(15)
    hold on;
    title('Green - diff 2')
end

figure(100)
hold on;
title('Green PPG')


% windowed operation
currDataWindow = rawData(1:numSamplesPerWindow , :);
currWindowStartSample = numSamplesPerWindow;
currTime = vecTime( 1:numSamplesPerWindow  );
if Ind_plotDebug
    [cwt_IR , f_IR] = cwt(rawData(:,2) , fs , 'VoicesPerOctave' , 24 , FrequencyLimits = Wn*(fs/2) );
    [cwt_G , f_G] = cwt(rawData(:,3) , fs , 'VoicesPerOctave' , 24 , FrequencyLimits = Wn*(fs/2) );
    figure
    surf( vecTime , f_IR , 20*log10(abs(cwt_IR)) , 'EdgeColor','none' );
    title('IR')
    figure
    surf( vecTime , f_G , 20*log10(abs(cwt_G)) , 'EdgeColor','none' );
    title('G')
end
vec_fundamentalFreq_estHR = zeros(1 , numTotalWindows);

for i = 1:numTotalWindows

    %% Step 1: Preprocessing
    % mean removal
    curr_DCcomponent = mean( currDataWindow , 1 );
    currDataWindow = currDataWindow - curr_DCcomponent;

    % filtering
    curr_RData = filter( b , a , currDataWindow(:,1) );
    curr_IRData = filter( b , a , currDataWindow(:,2) );
    curr_GData = filter( b , a , currDataWindow(:,3) );
    curr_filteredData = [curr_RData curr_IRData curr_GData];
    
    curr_ACcomponent =  max(curr_filteredData , [] , 1);
    curr_ACcomponent = curr_ACcomponent - min(curr_filteredData , [] , 1);

    %% Step 2: Motion Detection
    % MA component
    curr_MAcomponent = currDataWindow - curr_filteredData;

    % peak-to-noise Ratio for Green
    curr_filteredData_G_FFT = fftshift( fft( curr_filteredData(:,3) , NFFT , 1 ) );
    ind_Start = NFFT/2 + round(Wn(1)*NFFT/2);
    ind_End = NFFT/2 + round(Wn(2)*NFFT/2);
    peak = max( abs( curr_filteredData_G_FFT(ind_Start : ind_End) ) );
    sorted_freqDomain = sort( abs( curr_filteredData_G_FFT(ind_Start : ind_End) ) , 'ascend' );
    noise = mean( sorted_freqDomain( 1 : round( 0.2*(ind_End - ind_Start) ) ) );
    pToN_Ratio = peak / noise;

    % AC/DC ratio
    ACtoDC_Ratio = curr_ACcomponent ./ curr_DCcomponent;

    %% Step 3: CWT Motion removal
    % CWT
    [cwt_IR , f_IR] = cwt(curr_IRData , fs , 'VoicesPerOctave' , 24 , FrequencyLimits = Wn*(fs/2) );
    [cwt_G , f_G] = cwt(curr_GData , fs , 'VoicesPerOctave' , 24 , FrequencyLimits = Wn*(fs/2) );

    cwt_IR = cwt_IR ./ max( max( abs( cwt_IR ) ) );
    cwt_G = cwt_G ./ max( max( abs( cwt_G ) ) );
    
    % fundamental harmonic freq
    validInd = find( ((minFundamentalFreq/60) < f_G) & (f_G < (maxFundamentalFreq/60 + 0.25)) );

    [maxVal_G , maxInd_G] = max( sum( abs( cwt_G(validInd , :) ) , 2 ) );
    [maxVal_IR , maxInd_IR] = max( sum( abs( cwt_IR(validInd , :) ) , 2 ) );

    % isolation of harmonics and fundamental freq
    fundamentalFreq_estHR = f_G( validInd(maxInd_G) );
    vec_fundamentalFreq_estHR(i) = fundamentalFreq_estHR;
    numHarmonics = floor( f_G(1) / fundamentalFreq_estHR );
    harmonicInd = [];
    currMask = zeros( numel(f_G) , 1 );
    for ind_Harmonics = 1 : numHarmonics

        currInd = find( ( max( (fundamentalFreq_estHR*ind_Harmonics - stdDev_HR*3) , f_G(end) ) < f_G) & (f_G < min( (fundamentalFreq_estHR*ind_Harmonics + stdDev_HR*3) , f_G(1) ) ) );
        harmonicInd = [harmonicInd; currInd];

        currMask( currInd ) = 1;
        % currMask( currInd ) = normpdf(f_G(currInd) , fundamentalFreq_estHR*ind_Harmonics , stdDev_HR);        

    end

    % ICWT - reconstruction using harmonics
    cwt_G2 = cwt_G;
    cwt_G2 = cwt_G2 .* currMask;
    xrec_G = icwt(cwt_G2  , [] , f_G  , [f_G(end) f_G(1)] );
    curr_GData = icwt(cwt_G  , [] , f_G  , [f_G(end) f_G(1)] );
    if Ind_plotDebug
        [cwt_G_diff0 , f_G0] = cwt(xrec_G , fs , 'VoicesPerOctave' , 24 , FrequencyLimits = Wn*(fs/2) );
        figure(13)
        surf( currTime(1:end) , f_G0 , 20*log10(abs(cwt_G_diff0)) , 'EdgeColor','none' );
        hold off
    end

    % subtraction across remaining components
    % w_subtract = sum( abs( cwt_G(maxInd_IR , :) ) )  / sum( abs( cwt_IR(maxInd_IR , :) ) ) ;
    % cwt_subtract = cwt_G - w_subtract*cwt_IR;

    if Ind_plotDebug
        figure(10)
        surf( currTime , f_G , 20*log10(abs(cwt_IR)) , 'EdgeColor','none' );
        hold off
        figure(11)
        surf( currTime , f_G , 20*log10(abs(cwt_G)) , 'EdgeColor','none' );
        hold off
        % figure(12)
        % surf( currTime , f_G , 20*log10(abs(cwt_subtract)) , 'EdgeColor','none' );
        % hold off
    end

    % 

    %% Step 4: Temporal difference to enhance
    xrec_G_diff1 = diff( xrec_G );
    xrec_G_diff2 = diff( xrec_G_diff1 );
    
    if Ind_plotDebug
        [cwt_G_diff1 , f_G1] = cwt(xrec_G_diff1 , fs , 'VoicesPerOctave' , 24 , FrequencyLimits = Wn*(fs/2) );
        [cwt_G_diff2 , f_G2] = cwt(xrec_G_diff2 , fs , 'VoicesPerOctave' , 24 , FrequencyLimits = Wn*(fs/2) );
        figure(14)
        surf( currTime(2:end) , f_G1 , 20*log10(abs(cwt_G_diff1)) , 'EdgeColor','none' );
        hold off
        figure(15)
        surf( currTime(3:end) , f_G2 , 20*log10(abs(cwt_G_diff2)) , 'EdgeColor','none' );
        hold off
    end

    %% Plot
    figure(100);
    plot( currTime , curr_GData , 'b' );
    plot( currTime , xrec_G , 'r' );
    % plot( currTime(2:end) , xrec_G_diff1 , 'c' );
    % plot( currTime(3:end) , xrec_G_diff2 , 'm' );
    

    % next window samples
    if i < numTotalWindows
        currDataWindow = [ currDataWindow( 1 + numSamplesNewPerWindow : end , : ) ; rawData(currWindowStartSample+1:currWindowStartSample+numSamplesNewPerWindow , :) ];
        currTime = [ currTime( 1 + numSamplesNewPerWindow : end ) vecTime( currWindowStartSample+1:currWindowStartSample+numSamplesNewPerWindow ) ];
        currWindowStartSample = currWindowStartSample + numSamplesNewPerWindow;
    end

end
