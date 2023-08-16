% //////////////////////////////////////////////////////////////////////
% /////////////////////    Joukowskiy Airfoil      /////////////////////
% /////////////////////       Version 1.0       ////////////////////////
% Check README at https://github.com/DanielBelsky1999/JoukowskiyAirfoil
% //////////////////////////////////////////////////////////////////////


%            BE PATIENT!! The Code is quite slow! 
%            Please do NOT change the parameters too quickly. 

close all;
VERSION = "Release v1.0";

% Constants (Recommended settings! Other settings may cause unpredictable results)
C = 1;
U = 1;

% Circle Default Parameters
EPSILON = 0.14;
M = 0.18;
Gamma = 0;

% Zhukovskiy Transform (Recommended settings! Other settings may cause unpredictable results)
ZhukovskiyTransform = @(z, C) z+C/z;

UI_Init(C, EPSILON, M,U, Gamma, ZhukovskiyTransform, VERSION);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function UI_Init(C, EPSILON, M, U, Gamma, Transform, VERSION)
    
    TITLE = sprintf("Joukowskiy Transform %s", VERSION);
    
    fig_a = figure('Name',TITLE,'NumberTitle','off', "Resize", "off");
    fig_a.Units = "normalized";
    fig_a.Position = [0.3177, 0.1644, 0.4948, 0.7176];


    sgtitle("Joukowskiy Transform", "FontSize", 20);
    SubTitle1_Handle = uicontrol('Parent',figure(1),'Style','text',"Units","normalized", 'Position', [0.2368, 0.8548, 0.5263, 0.0484], "FontSize", 14);
    SubTitle2_Handle = uicontrol('Parent',figure(1),'Style','text',"Units","normalized", 'Position', [0.2368, 0.8065, 0.55, 0.0484], "FontSize", 12);

    axL = subplot(1,2,1);
    axL.Position = [0.13, 0.13, 0.3347, 0.8150];
    axR = subplot(1,2,2);
    axR.Position =  [0.5703, 0.13, 0.3347, 0.8150];
    axZoom = axes("Position", [0.79, 0.332, 0.115, 0.13]);
    axZoom.Visible = "off";
    plot_limits = [-3.5*C, 3.5*C, -3.5*C, 3.5*C];

    slider_eps = uicontrol('Parent',fig_a,'Style','slider','Units', 'normalized', 'Position', [0.1250, 0.2306, 0.3289, 0.0371],...
              'value',EPSILON, 'min',0, 'max',0.3, "SliderStep", [1/30, 0.1]);
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.0974, 0.2306, 0.0303, 0.0371],...
                'String','0');
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.4605, 0.2306, 0.0303, 0.0371],...
                    'String','0.3');
    EPSILON_Handle = uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.2289, 0.1855, 0.1316, 0.0371]);
    slider_M = uicontrol('Parent',fig_a,'Style','slider','Units', 'normalized', 'Position', [0.1250, 0.1452, 0.3289, 0.0371], ...
              'value',M, 'min',-0.5, 'max',0.5,  "SliderStep", [1/100, 0.1]);
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.0974, 0.1452, 0.0303, 0.0371],...
                'String','-0.5');
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.4605, 0.1452, 0.0303, 0.0371],...
                    'String','0.5');
    M_Handle = uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.2289, 0.1000, 0.1316, 0.0371]);
    slider_Gamma = uicontrol('Parent',fig_a,'Style','slider','Units', 'normalized', 'Position', [0.5737, 0.2323, 0.3289, 0.0371], ...
              'value',Gamma, 'min',-10, 'max',10,  "SliderStep", [1/100, 0.1]);
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.5382, 0.2306, 0.0303, 0.0371],...
                'String','-10');
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.9013, 0.2306, 0.0303, 0.0371],...
                    'String','10');
    GAMMA_Handle = uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.6763, 0.1855, 0.1316, 0.0371]); 
    check_box_StreamLines = uicontrol('style','checkbox','units','pixels',...
            'Units', 'normalized', 'Position', [0.5789, 0.1555, 0.0263, 0.0323]);
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.6053, 0.1524, 0.3947, 0.0371],...
                    'String','Stream Lines', "FontSize", 12,"HorizontalAlignment","left");
    check_box_Zoom = uicontrol('style','checkbox','units','pixels',...
            'Units', 'normalized', 'Position', [0.8, 0.1555, 0.0263, 0.0323]);
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.8264, 0.1524, 0.3947, 0.0371],...
                    'String','Zoom on T.E.', "FontSize", 12,"HorizontalAlignment","left");
    check_box_PressureMap = uicontrol('style','checkbox','units','pixels',...
            'Units', 'normalized', 'Position', [0.5789, 0.11, 0.0263, 0.0323], "Enable","off");
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.6053, 0.1070, 0.3947, 0.0371],...
                    'String','Pressure Map', "FontSize", 12,"HorizontalAlignment","left", "ForegroundColor", "#999999");
    check_box_Fill = uicontrol('style','checkbox','units','pixels',...
            'Units', 'normalized', 'Position', [0.8, 0.11, 0.0263, 0.0323]);
    uicontrol('Parent',fig_a,'Style','text','Units', 'normalized', 'Position', [0.8264, 0.1070, 0.3947, 0.0371],...
                    'String','Fill', "FontSize", 12,"HorizontalAlignment","left");
    button_Kutta = uicontrol('Parent',fig_a,'Style','pushbutton','Units', 'normalized', 'Position', [0.5526, 0.042, 0.3947, 0.0565],...
                    'String','Satisfy Kutta Condition', "FontSize", 12,"HorizontalAlignment","left");
    button_reset = uicontrol('Parent',fig_a,'Style','pushbutton','Units', 'normalized', 'Position', [0.2, 0.042, 0.2, 0.0565],...
                    'String','Reset', "FontSize", 12,"HorizontalAlignment","left");

    
    slider_eps.Callback = @(es,ed)    Complete_UI_Update(C, round(es.Value, 2), slider_M.Value, slider_Gamma.Value, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, check_box_StreamLines.Value, check_box_Zoom.Value, check_box_PressureMap.Value, check_box_Zoom.Value, axL, axR, axZoom, plot_limits);
    slider_M.Callback = @(es,ed)      Complete_UI_Update(C, slider_eps.Value, round(es.Value, 2), slider_Gamma.Value, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, check_box_StreamLines.Value,check_box_Zoom.Value, check_box_PressureMap.Value, check_box_Zoom.Value, axL, axR, axZoom, plot_limits);
    slider_Gamma.Callback = @(es,ed)  Complete_UI_Update(C, slider_eps.Value, slider_M.Value, round(es.Value, 2), U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, check_box_StreamLines.Value, check_box_Zoom.Value, check_box_PressureMap.Value, check_box_Zoom.Value,  axL, axR, axZoom, plot_limits);
    check_box_StreamLines.Callback = @(es, ed)    Complete_UI_Update(C, slider_eps.Value, slider_M.Value, slider_Gamma.Value, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, es.Value,check_box_Zoom.Value, check_box_PressureMap.Value, check_box_Fill.Value, axL, axR, axZoom, plot_limits);
    check_box_PressureMap.Callback = @(es, ed)    Complete_UI_Update(C, slider_eps.Value, slider_M.Value, slider_Gamma.Value, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, check_box_StreamLines.Value, es.Value,check_box_Zoom.Value, check_box_Fill.Value, axL, axR, axZoom, plot_limits);   
    check_box_Zoom.Callback = @(es, ed)    Complete_UI_Update(C, slider_eps.Value, slider_M.Value, slider_Gamma.Value, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, check_box_StreamLines.Value, es.Value, check_box_PressureMap.Value, check_box_Fill.Value, axL, axR, axZoom, plot_limits);   
    check_box_Fill.Callback = @(es, ed)    Complete_UI_Update(C, slider_eps.Value, slider_M.Value, slider_Gamma.Value, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, check_box_StreamLines.Value, check_box_Zoom.Value, check_box_PressureMap.Value, es.Value, axL, axR, axZoom, plot_limits);   
    button_Kutta.Callback = @(es, ed) SatisfyKuttaCondition(C, slider_eps.Value, slider_M.Value, U, slider_Gamma);
    button_reset.Callback = @(es, ed) Reset(C, EPSILON, M, Gamma, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, slider_eps, slider_M, slider_Gamma, check_box_StreamLines,check_box_Zoom, check_box_PressureMap, check_box_Fill, SubTitle1_Handle, SubTitle2_Handle, false,false, false,false, axL, axR, axZoom, plot_limits);

    Complete_UI_Update(C, EPSILON, M, Gamma, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, false,false, false,false, axL, axR, axZoom, plot_limits);
end

function Reset(C, EPSILON, M, Gamma, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, slider_eps, slider_M, slider_Gamma, check_box_StreamLines,check_box_Zoom, check_box_PressureMap, check_box_Fill, SubTitle1_Handle, SubTitle2_Handle, StreamLines_bool,Zoom_bool, PressureMap_bool,Fill_bool, axL, axR, axZoom, plot_limits)
    EPSILON_Handle.Value = EPSILON;
    M_Handle.Value = M;
    GAMMA_Handle.Value = Gamma;

    slider_eps.Value = EPSILON;
    slider_M.Value = M;
    slider_Gamma.Value = Gamma;
    check_box_StreamLines.Value = 0;
    check_box_PressureMap.Value = 0;
    check_box_Zoom.Value = 0;
    check_box_Fill.Value = 0;

    Complete_UI_Update(C, EPSILON, M, Gamma, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, StreamLines_bool,Zoom_bool, PressureMap_bool, Fill_bool, axL, axR, axZoom, plot_limits);
end

function Complete_UI_Update(C, EPSILON, M, Gamma, U, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle, StreamLines_bool,Zoom_bool, PressureMap_bool,Fill_bool, axL, axR, axZoom, plot_limits)
    [a, center, circle, airfoil] = Parameters_Init(C, EPSILON, M, Transform);
    UI_Update(C, EPSILON, M, Gamma, airfoil, Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle);
    Graphs_Update(a, C, EPSILON, M, Gamma, U, Transform, StreamLines_bool, Zoom_bool, PressureMap_bool, Fill_bool, center, circle, airfoil, axL, axR, axZoom, plot_limits);
end

function [a, center, circle, airfoil] = Parameters_Init(C, EPSILON, M, Transform)
        a = sqrt(M^2 + (C+EPSILON*C)^2);
        center = -EPSILON*C + M*1i;
        
        theta = 0:pi/50:2*pi;
        circle = center + a*exp(1i.*theta);
        airfoil = zeros(1, length(circle));
        for indx = 1:length(circle)
            airfoil(indx) = Transform(circle(indx), C);
        end
%         airfoil = [airfoil, airfoil(1)];
end

function UI_Update(C, EPSILON, M, Gamma, airfoil,  Transform, EPSILON_Handle, M_Handle, GAMMA_Handle, SubTitle1_Handle, SubTitle2_Handle)
        EPSILON_Handle.String = sprintf('Epsilon = %.2f', EPSILON);
        
        % Printed Value Correction
        if (M > -0.005)
            if (M < 0.005)
                M = 0;
            else
                M = abs(M);
            end
        end
        if (Gamma > -0.005)
            Gamma = abs(Gamma);
        end
        center_real = -C*EPSILON;
        if (center_real > -0.005)
            if (center_real < 0.005)
                center_real = 0;
            else
                center_real = abs(center_real);
            end
        end

        M_Handle.String = sprintf('M = %.2f', M);
        GAMMA_Handle.String = sprintf("Gamma = %.2f", Gamma);
        

        if(M>=0) 
            M_sign = "+"; 
        else 
            M_sign="-"; 
        end
        if (M == 0  &&  center_real == 0)
            SubTitle1_Handle.String = "Circle Center is at 0";
        elseif (M == 0)
            SubTitle1_Handle.String = sprintf("Circle Center is at %.2f", center_real);
        elseif (center_real == 0)
            SubTitle1_Handle.String = sprintf("Circle Center is at  %s %.2fi" ,M_sign, abs(M));
        else
            SubTitle1_Handle.String = sprintf("Circle Center is at %.2f %s %.2fi", center_real ,M_sign, abs(M));
        end

        chord_Length = Transform(C, C) - Transform(-EPSILON*C-C, C);
        SubTitle2_Handle.String = sprintf("Chord Length = %.2f, MaxThickness = %.2f", chord_Length, MaxThickness(C, airfoil));
end

function Graphs_Update(a, C, EPSILON, M, Gamma, U, Transform, StreamLines_bool, Zoom_bool, PressureMap_bool, Fill_bool, center, circle, airfoil, axL, axR, axZoom, plot_limits)

        black_axis_x = [plot_limits(1), plot_limits(2)];
        black_axis_y = [plot_limits(3), plot_limits(4)];     

        % SUBPLOT LEFT

        plot(axL, center,"r","Marker","+", "LineWidth", 0.6);
        hold(axL, "on");
        plot(axL, black_axis_x, zeros(1, length(black_axis_x)),"k");
        plot(axL, zeros(1, length(black_axis_y)), black_axis_y, "k");
    
        if (StreamLines_bool)
            stagnation_point_angle = asin(-Gamma/(4*pi*a*U));
            stagnation_point1 = a*exp(1i*stagnation_point_angle) - EPSILON*C + M*1i;
            stagnation_point2 = a*exp(1i*(pi-stagnation_point_angle)) - EPSILON*C + M*1i;
            plot(axL, [stagnation_point1,stagnation_point2],"LineStyle","none", "Color", "k", "Marker", "*");
            
            x = -4*C:0.01:4*C;
            y = -4*C:0.01:4*C;
            Psi_value_resolution = [-5:0.4:-1,-0.5,-0.15, 0.15,0.5,1:0.4:5];
            [X,Y] = meshgrid(x,y);
            Psi = U*sqrt((X).^2.+(Y).^2).*sin(atan2((Y),(X))).*(1.-((a^2)./((X).^2.+(Y).^2))) + (Gamma/(2*pi))*log(sqrt((X).^2.+(Y).^2));           
            for X_indx = 1:1:length(x)
                for Y_indx = 1:1:length(y)
                    if (sqrt((x(X_indx))^2 + (y(Y_indx))^2) <= a)
                        Psi(X_indx,Y_indx) = 0;
                    end
                end
            end
            x = x - C*EPSILON;
            y = y + M;
            [X,Y] = meshgrid(x,y);
            contour(axL, X,Y,Psi,Psi_value_resolution);
        end

%  Pressure Map (Button currently disabled)

        if (PressureMap_bool)
            x = -4*C:0.01:4*C;
            y = -4*C:0.01:4*C;
            [X, Y] = meshgrid(x,y);
            Velocity_at_XY = sqrt((U^2)*cos(1-(a./sqrt(X.^2+Y.^2)).^2).^2 + (U*sin(1+(a./sqrt(X.^2+Y.^2)).^2)-Gamma./(2*pi*sqrt(X.^2+Y.^2))).^2);
            Pressure = 0.5*U^2 - 0.5*Velocity_at_XY.^2;
            contourf(axL, X, Y, Pressure);
            colormap("turbo");
            fill(axR, real(airfoil), imag(airfoil), "k");
        end
    
        plot(axL, circle, "Color",[0.3,0,0.3]);
        	
        axL.DataAspectRatio = [1,1,1];
        axL.DataAspectRatioMode = "manual";
        axL.PlotBoxAspectRatioMode = "manual";
        
        grid(axL, "on");
        xlim(axL, [plot_limits(1), plot_limits(2)]);
        ylim(axL, [plot_limits(3), plot_limits(4)]);
        title(axL, "$\bf{\xi\ Plane}$", "Interpreter", "Latex");
        hold(axL, "off");
        
        % SUB PLOT RIGHT

        plot(axR, black_axis_x, zeros(1, length(black_axis_x)),"k");
        hold(axR, "on");
        plot(axR, zeros(1, length(black_axis_y)), black_axis_y, "k");
        
        if (StreamLines_bool)  
            plot(axR, [Transform(stagnation_point1, C),Transform(stagnation_point2, C)],"LineStyle","none", "Color", "k", "Marker", "*");

            z = zeros(length(x),length(y));
            for X_indx = 1:length(x)
                for Y_indx = 1:length(y)
                    z(Y_indx,X_indx) = Transform(x(X_indx)+1i*y(Y_indx), C);
                end
            end
            transformed_X = real(z);
            transformed_Y = imag(z);
            contour(axR, transformed_X, transformed_Y, Psi, Psi_value_resolution);

            % Zoom in
            if (Zoom_bool)
                zoom_in_Psi_value_resolution = [-0.3:0.01:-0.01, 0.01:0.01:0.3];
                contour(axZoom, transformed_X, transformed_Y, Psi, zoom_in_Psi_value_resolution);
                hold(axZoom, "on");
                trailing_edge = zeros(1, length(airfoil));
                for indx = 1:length(airfoil)
                    if (real(airfoil(indx)) > 1.7*C)
                        trailing_edge(indx) = airfoil(indx);
                    end
                end
                plot(axZoom, trailing_edge, "r", "LineWidth", 1.4);
                plot(axZoom, Transform(stagnation_point1, C), "Marker", "*", "Color","k");
                fill(axZoom, real(airfoil), imag(airfoil), "r");
                
                axZoom.XTick = [];
                axZoom.YTick = [];
    
                xlim(axZoom, [1.85*C, 2.05*C]);
                ylim(axZoom, [-0.1*C, 0.1*C]);
                hold(axZoom, "off");
            else
                axZoomOff(axZoom)
        
            end
        
        else
            axZoomOff(axZoom)

        end

%  Pressure Map (Button currently disabled) 
        if (PressureMap_bool)
            x = -4*C:0.01:4*C;
            y = -4*C:0.01:4*C;
            [X, Y] = meshgrid(x,y);
            z = zeros(length(x),length(y));
            for X_indx = 1:length(x)
                for Y_indx = 1:length(y)
                    z(Y_indx,X_indx) = Transform(x(X_indx)+1i*y(Y_indx), C);
                end
            end
            transformed_X = real(z);
            transformed_Y = imag(z);
            Velocity_at_XY = sqrt((U^2)*cos(1-(a./sqrt(X.^2+Y.^2)).^2).^2 + (U*sin(1+(a./sqrt(X.^2+Y.^2)).^2)-Gamma./(2*pi*sqrt(X.^2+Y.^2))).^2);
            Pressure = 0.5*(U^2 - Velocity_at_XY.^2);
            contourf(axR, transformed_X, transformed_Y, Pressure, [0.1, 0.2, 0.3]);
            colormap("turbo");
        end
        
        if (Fill_bool)
            plot(axR, airfoil, "k", "LineWidth", 0.6);
            fill(axR, real(airfoil), imag(airfoil), "r");
        else
            plot(axR, airfoil, "r", "LineWidth", 0.6);
        end

        axR.DataAspectRatio = [1,1,1];
        axR.DataAspectRatioMode = "manual";
        axR.PlotBoxAspectRatioMode = "manual";
        grid(axR, "on");
        xlim(axR, [plot_limits(1), plot_limits(2)]);
        ylim(axR, [plot_limits(3), plot_limits(4)]);
        title(axR, "$Z \bf{\ Plane}$", "Interpreter", "Latex");
        hold(axR, "off");
end

function axZoomOff(axZoom)
    if (axZoom.Visible == "on")
        axZoom.Visible = "off";
        for indx = 1:length(axZoom.Children)
            axZoom.Children(indx).Visible = "off";
        end
    end
end

function SatisfyKuttaCondition(C, EPSILON, M, U, slider_Gamma)

       a = sqrt(M^2 + (C+EPSILON*C)^2);
       M = round(M,2);
       theta = -asin(M/a);
       Calculated_Gamma = -4*pi*U*a*sin(theta);

   slider_Gamma.Value = Calculated_Gamma;
   es.Value = Calculated_Gamma;
   slider_Gamma.Callback(es, false);
end

function max_thickness = MaxThickness(C, airfoil)
    max_thickness = 0;
    %foreach point
    for indx = 1:length(airfoil)
        % check thickness closest points and calculate thickness
        for indx_2 = 1:length(airfoil)
            if (abs(real(airfoil(indx) - real(airfoil(indx_2)))) < 0.05*C)
                local_max = abs(imag(airfoil(indx)) - imag(airfoil(indx_2)));
                % compare to prev. max thickness
                if (local_max > max_thickness)
                    max_thickness = local_max;
                end
            end
        end 
    end
    % Set 0 if close to zero
    if max_thickness < 0.01*4*C
        max_thickness = 0;
    end
end