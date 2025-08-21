% =============================================================================
% MECHAIDE Matlab
% Version 1.0 [21-08-2025].
% -----------------------------------------------------------------------------
% File: mechaide_DXF.m v1.0 [21-08-2025].
% Description: Permite importar regiones de CAD .DXF
% =============================================================================
% Copyright (c) 2025 Nicolás Muñoz Guamán
% Email: nicolásjmunoz@gmail.com / nmunoz@mechaide.com
% Código original: www.tesis.mechaide.com
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% ============================================================================

classdef mechaide_DXF
    methods (Static)
        % Dibuja el contorno y traslada a 0,0
        function plot_data = plotNormal(DXF, paso_arco_inn, dividir_arco, salto_region)
            % Procesa datos DXF y genera puntos trasladados
            if nargin < 2
                paso_arco_inn = 2;
                dividir_arco = false;
                salto_region = true;
            end
            
            % Calcular minX y minY desde los datos originales
            minX = min([DXF(:,2); DXF(:,4)]);
            minY = min([DXF(:,3); DXF(:,5)]);
            
            % Aplicar traslado
            translatedData = mechaide_DXF.dataTranslate(DXF, minX, minY);
            
            % Generar puntos con la geometría trasladada
            plot_data = mechaide_DXF.plot(translatedData.array, paso_arco_inn, dividir_arco, salto_region);
        end
        
        % Traslada los datos DXF
        function translatedData = dataTranslate(DXF, cx, cy)
            array = DXF;
            array(:,2:7) = array(:,2:7) - [cx, cy, cx, cy, cx, cy];
            
            xmin = min([array(:,2); array(:,4)]);
            xmax = max([array(:,2); array(:,4)]);
            ymin = min([array(:,3); array(:,5)]);
            ymax = max([array(:,3); array(:,5)]);
            
            translatedData.array = array;
            translatedData.xmin = xmin;
            translatedData.xmax = xmax;
            translatedData.ymin = ymin;
            translatedData.ymax = ymax;
        end
        
        % Genera puntos a lo largo de las entidades DXF (líneas y arcos)
        function DXF_Aux = plot(DXF, paso_arco_inn, dividir_arco, salto_region)
            if nargin < 2
                paso_arco_inn = 2;
                dividir_arco = false;
                salto_region = true;
            end
            
            PI = pi;
            x = 0;
            DXF_Aux = {};
            idu = -1;
            id_old = 0;
            region_id_old = 0;
            
            for i = 1:size(DXF,1)
                paso_arco = paso_arco_inn;
                tipo_num = DXF(i,1);
                x0 = DXF(i,2);
                y0 = DXF(i,3);
                x1 = DXF(i,4);
                y1 = DXF(i,5);
                rx = DXF(i,6);
                ry = DXF(i,7);
                r = DXF(i,8);
                ars = DXF(i,9);
                are = DXF(i,10);
                id = DXF(i,11);
                region_id = DXF(i,12);
                giro = DXF(i,13);
                
                % Decodificar tipo_num a cadena
                switch tipo_num
                    case 1
                        tipo = 'LINE';
                    case 2
                        tipo = 'ARC';
                    case 3
                        tipo = 'CIRCLE';
                    otherwise
                        tipo = 'UNKNOWN';
                end
                
                if i < size(DXF,1)
                    region_id_fut = DXF(i+1,12);
                else
                    region_id_fut = region_id;
                end
                
                if strcmp(tipo, 'ARC') || strcmp(tipo, 'CIRCLE')
                    dife_angle = are - ars;
                    arco = abs(r * dife_angle);
                    
                    if dividir_arco
                        n_paso = floor(dife_angle * 13);
                    else
                        n_paso = floor(arco / paso_arco);
                    end
                    
                    if abs(dife_angle) > 1.25 * PI && n_paso < 5
                        n_paso = 6;
                    elseif abs(dife_angle) > 0.55 * PI && n_paso < 2
                        n_paso = 2;
                    end
                    
                    if n_paso == 0
                        paso_ang = dife_angle;
                        n_paso = 1;
                    else
                        paso_ang = dife_angle / n_paso;
                    end
                    
                    if giro == 1
                        for j = 0:n_paso
                            x0 = rx + r * cos(ars + j * paso_ang);
                            y0 = ry + r * sin(ars + j * paso_ang);
                            idu = idu + 1;
                            DXF_Aux{end+1} = struct('id', id, 'idu', idu, 'tipo', tipo, 'region_id', region_id, 'x', x0, 'y', y0);
                            x = x + 1;
                        end
                    else
                        for j = n_paso:-1:0
                            x1 = rx + r * cos(ars + j * paso_ang);
                            y1 = ry + r * sin(ars + j * paso_ang);
                            idu = idu + 1;
                            DXF_Aux{end+1} = struct('id', id, 'idu', idu, 'tipo', tipo, 'region_id', region_id, 'x', x1, 'y', y1);
                            x = x + 1;
                        end
                    end
                elseif strcmp(tipo, 'LINE')
                    idu = idu + 1;
                    DXF_Aux{end+1} = struct('id', id, 'idu', idu, 'tipo', tipo, 'region_id', region_id, 'x', x0, 'y', y0);
                    x = x + 1;
                    
                    if i == size(DXF,1) || region_id_fut ~= region_id
                        idu = idu + 1;
                        DXF_Aux{end+1} = struct('id', id, 'idu', idu, 'tipo', tipo, 'region_id', region_id, 'x', x1, 'y', y1);
                        x = x + 1;
                    end
                end
                
                region_id_old = region_id;
                id_old = id;
            end
        end
        
        % Procesa y plotea un archivo DXF con configuraciones simplificadas
        function plot_data = plotSimple(filename, paso_arco, plot, show_vertices)
            if nargin < 2
                paso_arco = 1; % Más puntos por arco por defecto
                show_vertices = false;
            end
            
            % Leer el archivo DXF
            rows = strsplit(fileread(filename), '\n');
            
            % Procesar
            DXF_raw = mechaide_DXF.read(rows);
            DXF_ordered = mechaide_DXF.order(DXF_raw);
            plot_data = mechaide_DXF.plotNormal(DXF_ordered.array, paso_arco, false, true);
            
            if (plot)
                % Separar puntos por región para evitar conectarlas
                x = zeros(length(plot_data), 1);
                y = zeros(length(plot_data), 1);
                regions = zeros(length(plot_data), 1);
                for i = 1:length(plot_data)
                    x(i) = plot_data{i}.x;
                    y(i) = plot_data{i}.y;
                    regions(i) = plot_data{i}.region_id;
                end
                
                % Plotear cada región por separado
                unique_regions = unique(regions);
                hold on;
                for r = 1:length(unique_regions)
                    region_mask = regions == unique_regions(r);
                    x_region = x(region_mask);
                    y_region = y(region_mask);
                    plot(x_region, y_region, 'k-', 'LineWidth', 1); % Líneas negras
                    if show_vertices
                        plot(x_region, y_region, 'k.', 'MarkerSize', 10); % Puntos gruesos para vértices
                    end
                end
                hold off;
                axis equal;
            end
        end
        
        % Procesa y plotea un archivo DXF
        function DXF = read(rows)
            DXF = {};
            Inicio = false;
            line_OK = false;
            x0_OK = false;
            y0_OK = false;
            x1_OK = false;
            y1_OK = false;
            arc_OK = false;
            rx_OK = false;
            ry_OK = false;
            r_OK = false;
            ars_OK = false;
            are_OK = false;
            circle_OK = false;
            rx_c_OK = false;
            ry_c_OK = false;
            r_c_OK = false;
            
            for x = 1:length(rows)
                line = strtrim(rows{x});
                if strcmp(line, 'ARC') || strcmp(line, 'LINE') || strcmp(line, 'CIRCLE')
                    Inicio = true;
                end
                if ~Inicio
                    continue;
                end
                if strcmp(line, 'ENDBLK')
                    return;
                end
                
                if x0_OK
                    x0 = str2double(line);
                    x0_OK = false;
                end
                if y0_OK
                    y0 = str2double(line);
                    y0_OK = false;
                end
                if x1_OK
                    x1 = str2double(line);
                    x1_OK = false;
                end
                if y1_OK
                    y1 = str2double(line);
                    y1_OK = false;
                    DXF = mechaide_DXF.xy(DXF, length(DXF), 'LINE', x0, y0, x1, y1);
                    line_OK = false;
                end
                
                if rx_OK
                    rx = str2double(line);
                    rx_OK = false;
                end
                if ry_OK
                    ry = str2double(line);
                    ry_OK = false;
                end
                if r_OK
                    r = str2double(line);
                    r_OK = false;
                end
                if ars_OK
                    ars = str2double(line) * pi / 180;
                    ars_OK = false;
                end
                if are_OK
                    are = str2double(line) * pi / 180;
                    are_OK = false;
                    if are < ars
                        are = are + 2 * pi;
                    end
                    x0 = rx + r * cos(ars);
                    y0 = ry + r * sin(ars);
                    x1 = rx + r * cos(are);
                    y1 = ry + r * sin(are);
                    DXF = mechaide_DXF.xy(DXF, length(DXF), 'ARC', x0, y0, x1, y1, rx, ry, r, ars, are);
                    arc_OK = false;
                end
                
                if rx_c_OK
                    rx_c = str2double(line);
                    rx_c_OK = false;
                end
                if ry_c_OK
                    ry_c = str2double(line);
                    ry_c_OK = false;
                end
                if r_c_OK
                    r_c = str2double(line);
                    r_c_OK = false;
                    x0 = rx_c + r_c * cos(0);
                    y0 = ry_c + r_c * sin(0);
                    DXF = mechaide_DXF.xy(DXF, length(DXF), 'CIRCLE', x0, y0, x0, y0, rx_c, ry_c, r_c);
                    circle_OK = false;
                end
                
                if strcmp(line, 'LINE')
                    line_OK = true;
                elseif strcmp(line, 'ARC')
                    arc_OK = true;
                elseif strcmp(line, 'CIRCLE')
                    circle_OK = true;
                end
                
                if line_OK && isnumeric(str2double(line))
                    code = str2double(line);
                    if code == 10
                        x0_OK = true;
                    elseif code == 20
                        y0_OK = true;
                    elseif code == 11
                        x1_OK = true;
                    elseif code == 21
                        y1_OK = true;
                    end
                elseif arc_OK && isnumeric(str2double(line))
                    code = str2double(line);
                    if code == 10
                        rx_OK = true;
                    elseif code == 20
                        ry_OK = true;
                    elseif code == 40
                        r_OK = true;
                    elseif code == 50
                        ars_OK = true;
                    elseif code == 51
                        are_OK = true;
                    end
                elseif circle_OK && isnumeric(str2double(line))
                    code = str2double(line);
                    if code == 10
                        rx_c_OK = true;
                    elseif code == 20
                        ry_c_OK = true;
                    elseif code == 40
                        r_c_OK = true;
                    end
                end
            end
        end
        
        % Ordena las entidades DXF
        function orderedData = order(DXF)
            PI = pi;
            first_ok = true;
            region_id = -1;
            x_first = 0;
            y_first = 0;
            tipo_old = ''; % Inicializar tipo_old
            x1_old = 0;    % Inicializar x1_old
            y1_old = 0;    % Inicializar y1_old
            
            % Preallocar array como matriz numérica
            array = zeros(length(DXF), 13); % 13 columnas: tipo, x0, y0, x1, y1, rx, ry, r, ars, are, id, region_id, giro
            object = {};
            
            for i = 1:length(DXF)
                tipo_str = DXF{i}{1};
                x0 = DXF{i}{2};
                y0 = DXF{i}{3};
                x1 = DXF{i}{4};
                y1 = DXF{i}{5};
                id = DXF{i}{6}; % id está en la posición 6 para LINE, 11 para ARC, 9 para CIRCLE
                
                % Inicializar valores predeterminados para ARC y CIRCLE
                rx = 0;
                ry = 0;
                r = 0;
                ars = 0;
                are = 0;
                
                % Asignar valores específicos para ARC y CIRCLE
                if strcmp(tipo_str, 'ARC')
                    rx = DXF{i}{6};
                    ry = DXF{i}{7};
                    r = DXF{i}{8};
                    ars = DXF{i}{9};
                    are = DXF{i}{10};
                    id = DXF{i}{11};
                elseif strcmp(tipo_str, 'CIRCLE')
                    rx = DXF{i}{6};
                    ry = DXF{i}{7};
                    r = DXF{i}{8};
                    id = DXF{i}{9};
                end
                
                % Codificar tipo como número: 1=LINE, 2=ARC, 3=CIRCLE
                switch tipo_str
                    case 'LINE'
                        tipo = 1;
                    case 'ARC'
                        tipo = 2;
                    case 'CIRCLE'
                        tipo = 3;
                    otherwise
                        tipo = 0; % Valor por defecto para casos inesperados
                end
                
                if i < length(DXF)
                    tipo_fut = DXF{i+1}{1};
                    x0_fut = DXF{i+1}{2};
                    y0_fut = DXF{i+1}{3};
                    x1_fut = DXF{i+1}{4};
                    y1_fut = DXF{i+1}{5};
                else
                    tipo_fut = '';
                    x0_fut = 0;
                    y0_fut = 0;
                    x1_fut = 0;
                    y1_fut = 0;
                end
                
                if strcmp(tipo_old, 'CIRCLE')
                    first_ok = true;
                end
                
                if strcmp(tipo_str, 'CIRCLE')
                    if i == length(DXF)
                        giro = 1;
                    else
                        giro = -1;
                    end
                    ars = 0;
                    are = 2 * PI;
                    first_ok = true;
                end
                
                if strcmp(tipo_str, 'LINE')
                    if first_ok
                        if (mechaide_DXF.round(x0_fut,2) == mechaide_DXF.round(x1,2) && mechaide_DXF.round(y0_fut,2) == mechaide_DXF.round(y1,2)) || ...
                                (mechaide_DXF.round(x1_fut,2) == mechaide_DXF.round(x1,2) && mechaide_DXF.round(y1_fut,2) == mechaide_DXF.round(y1,2))
                            giro = 1;
                        else
                            x0_temp = x0;
                            y0_temp = y0;
                            x0 = x1;
                            y0 = y1;
                            x1 = x0_temp;
                            y1 = y0_temp;
                            giro = -1;
                        end
                    elseif mechaide_DXF.round(x0,2) == mechaide_DXF.round(x1_old,2) && mechaide_DXF.round(y0,2) == mechaide_DXF.round(y1_old,2)
                        % No action needed
                    elseif mechaide_DXF.round(x1,2) == mechaide_DXF.round(x1_old,2) && mechaide_DXF.round(y1,2) == mechaide_DXF.round(y1_old,2)
                        x0_temp = x0;
                        y0_temp = y0;
                        x0 = x1;
                        y0 = y1;
                        x1 = x0_temp;
                        y1 = y0_temp;
                    end
                    giro = 0;
                end
                
                if strcmp(tipo_str, 'ARC')
                    if first_ok
                        if (mechaide_DXF.round(x0_fut,2) == mechaide_DXF.round(x1,2) && mechaide_DXF.round(y0_fut,2) == mechaide_DXF.round(y1,2)) || ...
                                (mechaide_DXF.round(x1_fut,2) == mechaide_DXF.round(x1,2) && mechaide_DXF.round(y1_fut,2) == mechaide_DXF.round(y1,2))
                            giro = 1;
                        else
                            x0_temp = x0;
                            y0_temp = y0;
                            x0 = x1;
                            y0 = y1;
                            x1 = x0_temp;
                            y1 = y0_temp;
                            giro = -1;
                        end
                    elseif mechaide_DXF.round(x0,2) == mechaide_DXF.round(x1_old,2) && mechaide_DXF.round(y0,2) == mechaide_DXF.round(y1_old,2)
                        giro = 1;
                    elseif mechaide_DXF.round(x1,2) == mechaide_DXF.round(x1_old,2) && mechaide_DXF.round(y1,2) == mechaide_DXF.round(y1_old,2)
                        x0_temp = x0;
                        y0_temp = y0;
                        x0 = x1;
                        y0 = y1;
                        x1 = x0_temp;
                        y1 = y0_temp;
                        giro = -1;
                    end
                end
                
                if strcmp(tipo_old, 'CIRCLE')
                    first_ok = true;
                end
                if first_ok
                    region_id = region_id + 1;
                    x_first = x0;
                    y_first = y0;
                    first_ok = false;
                elseif mechaide_DXF.round(x1,2) == mechaide_DXF.round(x_first,2) && mechaide_DXF.round(y1,2) == mechaide_DXF.round(y_first,2)
                    first_ok = true;
                end
                
                % Almacenar en matriz numérica
                array(i, :) = [tipo, x0, y0, x1, y1, rx, ry, r, ars, are, id, region_id, giro];
                object{end+1} = struct('tipo', tipo_str, 'x0', x0, 'y0', y0, 'x1', x1, 'y1', y1, ...
                    'rx', rx, 'ry', ry, 'r', r, 'ars', ars, 'are', are, 'id', id, ...
                    'region_id', region_id, 'giro', giro);
                
                tipo_old = tipo_str;
                x1_old = x1;
                y1_old = y1;
            end
            
            orderedData.array = array;
            orderedData.object = object;
        end
        
        % Redondea un valor a un número específico de decimales
        function val = rd(val, dec)
            if nargin < 2
                dec = 2;
            end
            val = round(val, dec);
        end
        
        % Redondea un valor (alias de rd)
        function val = round(val, dec)
            if nargin < 2
                dec = 2;
            end
            val = round(val, dec);
        end
        
        % Agrega una nueva entidad a la lista DXF
        function DXF = xy(DXF, j, tipo, varargin)
            new_row = {tipo, varargin{:}, j};
            DXF{end+1} = new_row;
        end
    end
end