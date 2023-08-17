 %% Remove curves

if curvasEliminadas ~= 0
    matrizInputElimin = matrizInput;
    matrizAbsoluteElimin = matrizAbsolute;
    %cont = 1;
    for i=1:length(matrizAbsolute(:,1))
        if (any(curvasEliminadas == i)) == 0
            matrizAbsoluteElimin(i,:) = matrizAbsolute(i,:);
           % cont = cont + 1;
        else
            if (any(curvasEliminadas == i)) == 1
                matrizAbsoluteElimin(i,:) = 0;
            end
        end
    end
    clear matrizAbsolute;
    matrizAbsolute = matrizAbsoluteElimin;
    clear matrizAbsoluteElimin;
end



%% Create individual matrices of each representation

for i=1:length(matrizAbsolute(:,1))
    matrizAbsoluteCortada(i,:) = matrizAbsolute(i,1:pointsAbs);
    SAbsolute = transpose(valoresS(1:pointsAbs,1));
    matrizHoltzerCortada(i,:) = matrizAbsolute(i,1:pointsHoltzer);
    SHoltzer = transpose(valoresS(1:pointsHoltzer,1));
    matrizKratkyCortada(i,:) = matrizAbsolute(i,1:pointsKratky);
    SKratky = transpose(valoresS(1:pointsKratky,1));
    SKratkyCuad = SKratky.*SKratky;
    matrizPorodCortada(i,:) = matrizAbsolute(i,1:pointsPorod);
    SPorod = transpose(valoresS(1:pointsPorod,1));
    SPorodCuat = SPorod.^4;
end

clear matrizAbsolute;
matrizAbsolute = matrizAbsoluteCortada;

for i=1:length(matrizHoltzerCortada(:,1))
    matrizHoltzer(i,:) = matrizHoltzerCortada(i,:) .* SHoltzer(1,:);
end
for i=1:length(matrizKratkyCortada(:,1))
    matrizKratky(i,:) = matrizKratkyCortada(i,:) .* SKratkyCuad(1,:);
end
for i=1:length(matrizPorodCortada(:,1))
    matrizPorod(i,:) = matrizPorodCortada(i,:) .* SPorodCuat(1,:);
end   
clear matrizAbsoluteCortada;
clear matrizHoltzerCortada;    
clear matrizKratkyCortada;
clear matrizPorodCortada;  
    
%% Scale matrices

% We do SVD of each matrix to obtain an adequate scale. We do with the
% clean and cutted data set (matrizKratky, matrizHoltzer y matrizPorod)

% SVD Holtzer
% for i=1:min(size(matrizHoltzer))
    [uHoltzer,sHoltzer,vHoltzer,xHoltzer, SigmaHoltzer] = pcarep(matrizHoltzer',min(size(matrizHoltzer)));
% end 
for i=1:length(sHoltzer)
    eigenvaluesHoltzer(i,1) = sHoltzer(i,i);
end

% SVD Kratky
% for i=1:min(size(matrizKratky))
    [uKrat,sKrat,vKrat,xKrat, SigmaKrat] = pcarep(matrizKratky',min(size(matrizKratky)));
% end
for i=1:length(sKrat)
    eigenvaluesKrat(i,1) = sKrat(i,i);
end

% SVD Porod
% for i=1:min(size(matrizPorod))
    [uPorod,sPorod,vPorod,xPorod, SigmaPorod] = pcarep(matrizPorod',min(size(matrizPorod)));
% end   
for i=1:length(sPorod)
    eigenvaluesPorod(i,1) = sPorod(i,i);
end

escaladoHoltzer = eigenvalues(1,1)/eigenvaluesHoltzer(1,1);
escaladoKratky = eigenvalues(1,1)/eigenvaluesKrat(1,1);
escaladoPorod = eigenvalues(1,1)/eigenvaluesPorod(1,1);


%% Create input matrix

if ((combinacionActual(1,1) == 1))
    matrizInput = matrizAbsolute;
    puntosMatrices = pointsAbs;
end
if ((combinacionActual(1,2) == 1))
    % Scale Holtzer
    matrizHoltzer = matrizHoltzer * escaladoHoltzer; 
    if ((combinacionActual(1,1) == 0))
        matrizInput = matrizHoltzer;
        puntosMatrices = pointsHoltzer;
    else
        matrizInput = [matrizInput matrizHoltzer];
        puntosMatrices = [puntosMatrices pointsHoltzer];
    end
end
if ((combinacionActual(1,3) == 1))
    % Scale Kratky
    matrizKratky = matrizKratky * escaladoKratky;
    matrizInput = [matrizInput matrizKratky];
    puntosMatrices = [puntosMatrices pointsKratky];
end

if ((combinacionActual(1,4) == 1))
    % Scale Porod
    matrizPorod = matrizPorod * escaladoPorod;
    matrizInput = [matrizInput matrizPorod];
    puntosMatrices = [puntosMatrices pointsPorod];
end

if InitEstmethod == 2
%% Create initial estimations dataset 
firstcurve = horzcat(valoresS, Intensities(:,1), Errors(:,1));
for i = 1:length(Intensities(1,:))
curvetocompare = horzcat(valoresS,Intensities(:,i));
chiforcomp1(i,1) = i;
chiforcomp1(i,2) = compare2curves(firstcurve, curvetocompare);
end
[~,idx] = sort(chiforcomp1(:,2));
sortedchi1values = chiforcomp1(idx,:);
secondcurve = horzcat(valoresS, Intensities(:,TwoEstimates(2)), Errors(:,TwoEstimates(2)));
for i = 1:length(Intensities(1,:))
curvetocompare = horzcat(valoresS, Intensities(:,i));
chiforcomp2(i,1) = i;
chiforcomp2(i,2) = compare2curves(secondcurve, curvetocompare);
end
[~,idx] = sort(chiforcomp2(:,2));
sortedchi2values = chiforcomp2(idx,:);

for i = 1:length(sortedchi1values)
curvenumber = sortedchi1values(i,1);
[~,idx]=min(abs(curvenumber-sortedchi2values(:,1)));
chirankavglist(i,1) = curvenumber;
chirankavglist(i,2) = (i+idx)/2;
end
[~,idx] = sort(chirankavglist(:,1));
sortedchirankavg = chirankavglist(idx,:);
n = 1;
chiforcompavg = (chiforcomp1 + chiforcomp2)/2;
% [M,I] = max (chiforcompavg);
% chicompcomb(:,1) = chiforcomp1;
% chicompcomb(:,2) = chiforcomp2;
% chistddev(:,1) = std(chicompcomb,0,2);
[M,I] = max (sortedchirankavg(:,2));
ThirdCurveNumber = I;

PureComponents = horzcat(1,TwoEstimates(2),I);
end
if InitEstmethod == 1
    PureComponents = curvasEstInic;
end
if InitEstmethod == 3
    PureComponents = Estimates;
end
curvasEstInicAbs = PureComponents;
curvasEstInicHoltzer = PureComponents;
curvasEstInicKratky = PureComponents;
curvasEstInicPorod = PureComponents;

if InitEstmethod == 3
    PureComponents = Estimates;
end
curvasEstInicAbs = PureComponents;
curvasEstInicHoltzer = PureComponents;
curvasEstInicKratky = PureComponents;
curvasEstInicPorod = PureComponents;



for i=1:numeroEspecies
    if ((combinacionActual(1,1) == 1))
        estimInicialesAbs(i,:) = matrizAbsolute((curvasEstInicAbs(1,i)),:);
    end
    if ((combinacionActual(1,2) == 1))
        estimInicialesHoltzer(i,:) = matrizHoltzer((curvasEstInicHoltzer(1,i)),:);
    end 
    if ((combinacionActual(1,3) == 1))
        estimInicialesKratky(i,:) = matrizKratky((curvasEstInicKratky(1,i)),:);
    end   
    if ((combinacionActual(1,4) == 1))
        estimInicialesPorod(i,:) = matrizPorod((curvasEstInicPorod(1,i)),:);
    end   
end
if ((combinacionActual(1,1) == 1))
    estimInicialesInput = estimInicialesAbs;
end
if ((combinacionActual(1,2) == 1))
    estimInicialesInput = [estimInicialesInput estimInicialesHoltzer];
end 
if ((combinacionActual(1,3) == 1))
    estimInicialesInput = [estimInicialesInput estimInicialesKratky];
end
if ((combinacionActual(1,4) == 1))
    estimInicialesInput = [estimInicialesInput estimInicialesPorod];
end


%% Create Equality constrain for spectra (if exist)
clear ssel;
if curvasFijadas == 1
    
    for ll=1:numeroEspecies
        if curvasFijadas2(1,ll) == 0
            ssel(ll,1:length(matrizInput)) = NaN;
        end
        if curvasFijadas2(1,ll) == 1
            ssel(ll,:) = matrizInput(curvasEstInic(1,ll),:);
        end
    end
    

else
    ssel = 0;
end




%% Define number of matrices 

numMatrices = 0;

for i=1:length(combinacionActual(1,:))
    if combinacionActual(1,i) == 1
        numMatrices = numMatrices + 1;
    end
end



