%% setup paths


clear
restoredefaultpath
addpath(genpath('submodules'));
addpath('~/git/gisogrimm/tascar/scripts')
check_for_tascar()



load Thiemann2019.mat
eval_az_deg = az_deg;
load thiemanngt.mat


[b,a]=v_stdspectrum(2,'z',44100);

%% Find approximation
%plugins =["finalresponsensp.mat","finalresponsevbap.mat","finalresponsehoa.mat"];
vspknos = [-1, -1,0];%-1 for using kernel no. 0 for no vspks, any other no. is a constant number of vspks

figure(1)
figure(2)
figure(3)
for n = 1:length(plugins)
    plugin = plugins(n)
    load(plugin);
    inperr = zeros(1,length(spkno));
    recerr =zeros(1,length(spkno));
    worstinperr = zeros(1,length(spkno));
    worstdinperr = zeros(1,length(spkno));
    powerr = zeros(1,length(spkno));
    for i = 1:length(spkno)
        ttruth = truth;
        if (~contains(plugins(n),"vbap") && ~contains(plugins(n),"nsp") && ~contains(plugins(n),"vbip"))
            disp("shifting by 1");
            ttruth = truth(2:512,:,:);
            rresponse = squeeze(response(i,1:511,:,:));
        else
            rresponse = squeeze(response(i,:,:,:));
        end

        %now separate directions where vspk used
        vspkno =0;
        if vspknos(n) == -1
            vspkno = spkno(i);
        elseif vspknos(n) > 0
            vspkno = vspknos(n);
        end
        vspkdirs = [];
        for spk = 1:vspkno
            dir = (spk-1)*720/vspkno;
            if dir == round(dir)
                vspkdirs = [vspkdirs, round(dir)+1];
            end
        end
        
        rec = rresponse(:,:,vspkdirs);
        rect =ttruth(:,:,vspkdirs);
        inp = rresponse;
        inp(:,:,vspkdirs) = [];
        inpt = ttruth;
        inpt(:,:,vspkdirs) = [];
        if ~isempty(rec)
            recerr(i) = sum(((rec-rect).^2),'all')/numel(eval_az_deg);
        end
        err = (filter(b,a,inp-inpt)).^2;
        
        %errarr(i) = sum(err, 'all')/numel(err);
        tnorm = sqrt(sum(filter(b,a,inpt).^2,1));
        terr = sqrt(sum(err,1))./tnorm;
        inperr(i) = sum(terr,'all')/numel(terr);
        worstdinperr(i)=max(terr,[],'all');

    end

    
    figure(1)
    plot(spkno,mag2db(inperr))
    hold on
    
    if vspknos(n)~= 0
        figure(2)
        plot(spkno,recerr)
        hold on
    end

    figure(3)
    plot(spkno,mag2db(worstdinperr))
    hold on
    
end
figure(1)
legend(legtext)
xlabel("Number of kernels")
ylabel("Mean error in impulse response (dB)")
figure(2)
legend(legtext)
xlabel("Number of kernels")
ylabel("Mean reconstruction error in impulse response")
figure(3)
legend(legtext)
xlabel("Number of kernels")
ylabel("Worst error in impulse response (dB)")





