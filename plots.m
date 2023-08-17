cc = jet(length(curvesFinal));

if seleccion == 1;

    figure();
    subplot(2,2,1);
    for i=1:curvasTotales
        plot(valoresS(1:pointsAbs,1), Intensities(1:pointsAbs,i),'Color',cc(i,:), 'LineWidth',2);
        hold on;
        ylabel ('I(q)', 'FontSize',15)
        xlabel ('q', 'FontSize',15)
    end
    title('Absolute scale');

    subplot(2,2,2);
    for i=1:curvasTotales
        plot(valoresS(1:pointsHoltzer,1), Intensities(1:pointsHoltzer,i).*valoresS(1:pointsHoltzer,1),'Color',cc(i,:), 'LineWidth',2);
        hold on;
        ylabel ('I(q)*q', 'FontSize',15)
        xlabel ('q', 'FontSize',15)
    end
    title('Holtzer plot');

    subplot(2,2,3);
    valoresSsquare = valoresS.*valoresS;
        for i=1:curvasTotales
        plot(valoresS(1:pointsKratky,1), Intensities(1:pointsKratky,i).*valoresSsquare(1:pointsKratky,1),'Color',cc(i,:), 'LineWidth',2);
        hold on;
        ylabel ('I(q)*q^2', 'FontSize',15)
        xlabel ('q', 'FontSize',15)
    end
    title('Kratky plot');

    subplot(2,2,4);
    valoresScuad = valoresS.^4;
    for i=1:curvasTotales
        plot(valoresS(1:pointsPorod,1), Intensities(1:pointsPorod,i).*valoresScuad(1:pointsPorod,1),'Color',cc(i,:), 'LineWidth',2);
        hold on;
        ylabel ('I(q)*q^4', 'FontSize',15)
        xlabel ('q', 'FontSize',15)
    end
    title('Porod plot');

end

if seleccion == 2
    
    
    
    figure();
    for i=1:length(curvesFinal)
        semilogy(valoresS(:,1),Intensities(:,i),'Color',cc(i,:), 'LineWidth',2);
        ylabel ('Intensity', 'FontSize',15)
        xlabel ('q', 'FontSize',15)
        hold on;
    end
    title('Data set in semilogarithmic scale');


end


if seleccion == 3
    
    
    
    figure();
    for i=1:length(extraData)
        plot(intensidad(:,i),'Color',cc(i,:));
        hold on;
    end
    title('Additional Data Set');


end



