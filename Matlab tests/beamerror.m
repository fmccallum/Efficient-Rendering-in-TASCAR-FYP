%% setup paths


clear
restoredefaultpath
addpath(genpath('submodules'));
addpath('~/git/gisogrimm/tascar/scripts')
check_for_tascar()

addpath('submodules/sap-elobes-microphone-arrays');
addpath('submodules/sap-elobes-utilities');
addpath('submodules/sap-voicebox/voicebox');
addpath('lib');


spk_file = 'temp.spk';

load Thiemann2019.mat
load thiemanngt.mat
eval_az_deg = az_deg;
%% beamformer design
% ensure empty struct
spec = struct();

% sample rate
spec.fs = fs;            

% cartesian unit vector
spec.look_unit_vec = [1 0 0];  

% microphone array
nSensors = 5;
spacing = 0.03;
%spec.ema_fcn = @()FreeField_20180202_01_ULA_configurable(nSensors,spacing);
spec.ema_fcn =@Sampled_20220225_Thiemann2019_bte_horiz;

filtCoefs = design_beamformer(spec);

ema = spec.ema_fcn();
%thiemann ema.setSoundSpeed(v_soundspeed);
ema.prepareData(spec.fs);
[target_az,target_inc,~] = mycart2sph(spec.look_unit_vec);

%%These bits
h = ema.getImpulseResponseForSrc(target_az,target_inc);

%per channel filtered response
filtered = fftfilt(filtCoefs,[h;zeros(size(filtCoefs,1)-1,ema.nSensors)]);
%Beamformer output
beamout = sum(filtered,2);


vspk = struct();
vspk.ema_fcn = spec.ema_fcn;
vspk.fs = spec.fs;

type = 'nsp'


%% Get true response

N = length(eval_az_deg);
energypattern = zeros(802,N);
trueenergypattern = zeros(802,N);
[b,a]=v_stdspectrum(2,'z',spec.fs);


for idoa = N:-1:1
    h = truth(:,:,idoa);

    %per channel filtered response
    filtered = fftfilt(filtCoefs,[h;zeros(size(filtCoefs,1)-1,ema.nSensors)]);
    %Beamformer output
    beamout = sum(filtered,2);
    trueenergypattern(:,idoa) = filter(b,a,beamout);
end


figure;
polarplot(deg2rad(eval_az_deg),trueenergypattern')
title("RMS of A weighted filtered response");



%% Find approximation
plugins =["finalresponsensp.mat","finalresponsevbap.mat","finalresponsehoa.mat"];
plugins=["vbapresponse.mat","vbipresponse.mat"]


figure;
for n = 1:length(plugins)
    plugin = plugins(n)
    load(plugin);
    energypattern = zeros(802,N);
    errarr = zeros(1,length(spkno));
    for i = 1:length(spkno)
        spkno(i)
        y = squeeze(response(i,:,:,:));
        for idoa = N:-1:1
            h = y(:,:,idoa);

            %per channel filtered response
            filtered = fftfilt(filtCoefs,[h;zeros(size(filtCoefs,1)-1,ema.nSensors)]);
            %Beamformer output
            beamout = sum(filtered,2);
            energypattern(:,idoa) = filter(b,a,beamout);
        end
        if (~contains(plugins(n),"vbap") && ~contains(plugins(n),"nsp") && ~contains(plugins(n),"vbip"))
            err = abs(trueenergypattern(2:end)-energypattern(1:end-1)).^2;
        else
            err = abs(trueenergypattern-energypattern).^2;
        end    
        terr =sum(err, 1);
        errarr(i) = sqrt(sum(terr, 'all'))/numel(terr);
    end
    semilogx(spkno,mag2db(errarr))
    hold on;

end

legend(["NSPK","VBAP","SH  "],'Location','east')
xlabel("Number of kernels")
ylabel("Worst case error in beam pattern (dB)")





