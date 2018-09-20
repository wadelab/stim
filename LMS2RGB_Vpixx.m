function [stimRGB, maxcontrast] = LMS2RGB_Vpixx (stimLMS, fundamentals, resampledSpectra)

% Attempt to get RGB triplet for a specified LMS excitation.
% Using spectra obtained on the Viewpixx and through the goggles, 2/4/16
% or, spectra for the PROpixx obtained on .... 
% That have been resampled to fit the same set of wavelengths as the cone fundamentals.
% The result will depend on the spectra that are loaded (ie, Viewpixx or PROpixx).
% The Viewpixx spectra are loaded by default.
% Make sure the resampled spectra are in the proper format...
%
% Based on cone2RGB & findMaxConeScale from the vistadisp tools suite.
% Also return the max possible contrast (RTM, 13/5/16). Note that this value is the 
% same for negative & positive polarities.
% R Maloney, April 2016

%------------------------------------------

%From cone2RGB....

% ARGUMENTS

%  stimLMS:  .dir    is the color direction of the contrast stimulus.
%             By convention,
%             .dir is a 3-value vector with a maximum value of 1.0, for L, M, S.
%            .scale  is the scale factor (ie amount of desired cone contrast).
%            When the stimLMS.dir is cone isolating, the
%            scale factor is the same as contrast.  The
%            definition of a single contrast value is problematic for other
%            directions. 
% fundamentals: optional: the cone fundamentals. So far we have been using Stockman & Sharpe 2000
%                2 deg cone fundamentals sample every 1 nm.
%               If not entered, they are loaded here.
%               But we could also use the 10 deg cone fundamentals; it depends on what is loaded.
% spectra:     The RESAMPLED monitor/projector spectra. 
%              For the Viewpixx these were obtained 2/4/16.
%              For the PROpixx, they were obtained ... 
%              They come from the Ocean Optics Spectra Suite software using the Jaz spectrometer.
%              They have been resampled to integer increases in nm, so they can be associated with 
%              The same wavelengths in the cone fundamentals (also sampled at 1 nm integer steps).
%              spectra takes the form of a 1 *3 cell array.
%              Cells 1:2 are each of size 2048 4 & contain the spectra for R,G,B & K guns (at max luminance).
%              The two cells are for the left and right eye data, respectively.
%              Cell {3} contains the wavelengths for all the spectra.

% RETURNS
%
% stimRGB: 
%          .dir is set to rgb direction corresponding to the desired lms direction. 
%          .scale is what you scale the RGB .dir by, kind of like contrast?
%          .maxScale give you the max. scale factor possible with your monitor gamut.
%          an error will occur if your requested scale factor is larger than the max scale.
%          Once we have stimRGB, to get your final RGB triplet for drawing your stimulus,
%          evaluate stimRGB.dir * stimRGB.scale + 0.5 (the background)

%set the desired LMS coordinates, for example,:
%stimLMS.dir = [0 0 1]; % for S-cone
%stimLMS.scale = 0.6; %amount of cone contrast.

%Set the background 
%(optional) .dir and  .scale define the mean RGB of background,
% so that backRGB.dir*backRGB.scale is a vector of linear rgb values.
% We will just use the default here:
%disp('LMS2RGB: Using default background of [0.5 0.5 0.5]')
backRGB.dir = [1 1 1]';
backRGB.scale = 0.5;

% --------------------------------------------------------------- %
% load the stockman/sharpe 2 deg cone fundamentals if not parsed
% --------------------------------------------------------------- %
if nargin < 2 || isempty(fundamentals)
    load('\\PSHome\Home\rm1380\My Documents\Colour\StockmanSharpe_2deg_cone_fundamentals_1nm.mat')
end
sensors = fundamentals(:,2:4); %cut the first column, which gives wavelength

% --------------------------------------------------------------- %
%                   Load monitor spectra
% --------------------------------------------------------------- %

if nargin < 3 || isempty (resampledSpectra)
    % Load the Viewpixx spectra, acquired 2/4/16
    % Note that these include both raw and resampled spectra.
    % We want the spectra that have been resampled to a range increasing by 1 integer nm (as do the cone fundamentals)
    load('\\PSHome\Home\rm1380\My Documents\Calibration\Viewpixx_Processed_cal_data_2_4_2016.mat')
end
desiredwl = resampledSpectra{3}; %pull out the wavelengths for the spectra.

%We want to restrict the spectra to the same range of wavelengths as the cone fundamentals
coneRange = [fundamentals(1,1), fundamentals(end,1)];
spectRange = [find(desiredwl==coneRange(1)), find(desiredwl==coneRange(2))];

% We will take the mean of the left and right eye spectra, because they are so similar.
spectra = mean([resampledSpectra{1}(:,1),resampledSpectra{2}(:,1)], 2); %Mean of L & R eye Red spectra
spectra(:,2) = mean([resampledSpectra{1}(:,2),resampledSpectra{2}(:,2)], 2); %Mean of L & R eye Green spectra
spectra(:,3) = mean([resampledSpectra{1}(:,3),resampledSpectra{2}(:,3)], 2); %Mean of L & R eye Blue spectra

%Now restrict those spectra to the same range as the cone fundamentals
spectra = spectra(spectRange(1):spectRange(2),:);

%The spectra and fundamentals (sensors) matrices now are the same size.
% Plot the spectra & the fundamentals out if you wish.
% note that both plots now use the same x axis values.
% figure
% subplot(1,2,1)
% plot(fundamentals(:,1), fundamentals(:,2), 'r-')
% hold on
% plot(fundamentals(:,1), fundamentals(:,3), 'g-')
% plot(fundamentals(:,1), fundamentals(:,4), 'b-')
% ylabel('Relative sensitivity')
% xlabel('Wavelength (nm)')
% subplot(1,2,2)
% plot(fundamentals(:,1), spectra(:,1), 'r-') 
% hold on
% plot(fundamentals(:,1), spectra(:,2), 'g-') 
% plot(fundamentals(:,1), spectra(:,3), 'b-') 
% ylabel('Absolute irradiance (uW/cm^2/nm)')
% xlabel('Wavelength (nm)')

% Now, multiply the 2 matrices to get the rgb2lms conversion matrix:
rgb2lms = sensors'*spectra;

% At this point, cone2RGB.m calls findMaxConeScale.m, so we will pick things up there:
%---------------------------------------

% [stimLMS, stimRGB] = findMaxConeScale(rgb2lms,stimLMS,backRGB);

%--------------------------------------
% From findMaxConeScale:

%  dirFlag:  if true, we assume the stimulus is a positive modulation from
%              the background. If 0 or unset, we assume the stimulus
%              modulates symmetrically in two directions from the
%              background. This affects how we calculate the maximum scale
%              along a color direction.
%if notDefined('dirFlag'), dirFlag = 0; end
dirFlag = 0;

%Compute the inverse matrix of rgb2lms:
lms2rgb = inv(rgb2lms);

% Check whether the background RGB values are within the unit cube
meanRGB = backRGB.dir(:) * backRGB.scale; %mean RGB: a mean luminance screen
err = checkRange(meanRGB,[0 0 0]',[1 1 1]');
if err ~= 0,  error('meanRGB out of range'); end

%  Determine the background LMS direction   
lmsBack = rgb2lms*(meanRGB);

%  Scale stimulus LMS by the background LMS
%  We do this because the stimLMS dir is a contrast. So if we want a 0.1
%  contrast in the L and a 0.2 contrast in the M, then the amount of L and
%  M we need would be .1*lmsBack(1) and .2*lmsBack(2).
scaledStimLMS = stimLMS.dir(:) .* lmsBack;

%  Determine the stimulus RGB direction that will create the desired LMS
%  changes.  We should probably not create lms2rgb as above, but use the
%  \ operator instead.
stimRGB.dir = lms2rgb*scaledStimLMS; %#ok<MINV>
stimRGB.dir = stimRGB.dir/max(abs(stimRGB.dir));

% We want to find the largest scale factor such that the
% background plus stimulus fall on the edges of the unit cube.
% We begin with the zero sides of the unit cube, 
% 
%      zsFactor*(stimRGB.dir) + meanRGB = 0
% 
% Solving this equation for zsFactor, we obtain
%
sFactor = -(meanRGB) ./ stimRGB.dir;

%  The smallest scale factor that bumps into this side is
if dirFlag
    zsFactor = min(abs(sFactor(sFactor > 0))); % if we check only positive modulations from background
else
    zsFactor = min(abs(sFactor)); % if we assume symmetric, bidirectional modulation from background
end          

% Now find the sFactor that limits us on the 1 side of the unit RGB cube.
% 
%       usFactor*stimRGB.dir + meanRGB = 1
%   
sFactor = (ones(3,1) - meanRGB) ./ stimRGB.dir;
 
if dirFlag
    usFactor = min(abs(sFactor(sFactor>0))); % if we check only positive modulations from background
else
    usFactor = min(abs(sFactor));    % if we assume symmetric, bidirectional modulation from background
end


%  Return the smaller of these two factors
%  
stimRGB.maxScale = min(zsFactor,usFactor);

% Next, convert these values into LMS contrast terms.
% 
% General discussion:
% 
%  For each scale factor applied to the stimulus, there is a
%  corresponding contrast.  But, this must be computed using both
%  the stimLMS and the backLMS.  So, contrast and stimLMS.scale
%  are not uniquely linked, but they depend on the background.
% 
%  When stimRGB.scale is less than stimRGB.maxScale, we are sure that we
%  are within the unit cube on this background.  What is the
%  highest scale level we can obtain for the various cone classes
%  at this edge? 
% 
% Compute the LMS coordinates of the [stimulus plus background] and
% the background alone.  Use these to compute the max scale
% factor we can use in the LMS direction.  This is the maximum
% contrast when we are in a cone isolating direction.
%  
lmsStimPlusBack = ...
    rgb2lms*(stimRGB.maxScale*stimRGB.dir + backRGB.dir*backRGB.scale);
lmsContrast = (lmsStimPlusBack  - lmsBack) ./ lmsBack;
stimLMS.maxScale = max(abs(lmsContrast));

maxcontrast = norm(lmsContrast); %this value is provided by findMaxConeScale but not actually asked for by cone2RGB.

% *** the end of the stuff from findMaxConeScale.m ***

%------------------------------------
% Now the stuff from coneRGB again, where it left off, now that we have the variables stimLMS & stimRGB:
for ii=1:length(stimLMS.scale)
	if (stimLMS.scale(ii) > stimLMS.maxScale)
		if (stimLMS.scale(ii)-stimLMS.maxScale < 0.001)
			stimLMS.scale(ii) = stimLMS.maxScale;
		else
	      	%error('Requested contrast ( %.3f) exceeds maximum (%.3f)\n', stimLMS.scale(ii),stimLMS.maxScale);
            % Instead of giving an error, we will just make the stimLMS.scale = stimLMS.maxScale
            % This ultimately makes the initial contrast level irrelevant.
            % If you want max contrast, set it to 1. It will set it to the max. possible, & give a message it has done so.
            % This actually makes varying the contrast rather easy.
            % We had to change this bit because we couldn't have it crashing every time the values changed in the min flicker isoluminance task.
            % RM - 22/4/16
            fprintf('Requested contrast ( %.3f) set to maximum (%.3f)\n', stimLMS.scale(ii),stimLMS.maxScale);
            stimLMS.scale(ii) = stimLMS.maxScale;
		end  
	end
    % When stimRGB.scale equals stimRGB.maxScale, 
    % 
    %      stimLMS.scale = stimLMS.maxScale
    % 
    % Everything is linear, so to obtain 
    % 
    %    stimLMS.scale = stimLMS.maxScale * (stimRGB.scale/stimRGB.maxScale)
    % 
    % To solve for the stimRGB.scale that yields a stimLMS.scale,
    % we invert
    stimRGB.scale(ii) = (stimLMS.scale(ii)/stimLMS.maxScale)*stimRGB.maxScale;
end


% Debugging
% Compute the stimulus contrast for the various conditions, as
% a check.
% (these values also appear to be unused).
lmsBack = rgb2lms*(backRGB.dir*backRGB.scale);
for ii=1:length(stimRGB.scale)
  lmsStimPlusBack = rgb2lms*(stimRGB.scale(ii)*stimRGB.dir) + lmsBack;
  (lmsStimPlusBack - lmsBack) ./ lmsBack;
end


%--------------------------------------------
% In case you'd like to plot the spectra, each eye:
% figure
% %Left eye:
% subplot(1,2,1)
% plot(desiredwl(spectRange(1):spectRange(2)), resampledSpectra{1}(spectRange(1):spectRange(2),1), 'r')
% hold on
% plot(desiredwl(spectRange(1):spectRange(2)), resampledSpectra{1}(spectRange(1):spectRange(2),2), 'g')
% plot(desiredwl(spectRange(1):spectRange(2)), resampledSpectra{1}(spectRange(1):spectRange(2),3), 'b')
% xlabel('wavelength (nm)')
% ylabel('Absolute Irradiance (uW/cm^2/nm)')
% title('Left eye, through goggles')
% % Right eye:
% subplot(1,2,2)
% plot(desiredwl(spectRange(1):spectRange(2)), resampledSpectra{2}(spectRange(1):spectRange(2),1), 'r')
% hold on
% plot(desiredwl(spectRange(1):spectRange(2)), resampledSpectra{2}(spectRange(1):spectRange(2),2), 'g')
% plot(desiredwl(spectRange(1):spectRange(2)), resampledSpectra{2}(spectRange(1):spectRange(2),3), 'b')
% xlabel('wavelength (nm)')
% ylabel('Absolute Irradiance (uW/cm^2/nm)')
% title('Right eye, through goggles')
% suptitle('PROpixx spectra, 21/4/16')