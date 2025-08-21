% =============================================================================
% MECHAIDE Matlab
% Version 1.0 [21-08-2025].
% -----------------------------------------------------------------------------
% File: mechaide_MESH.m v1.0 [21-08-2025].
% Description: Genera malla DPOLYLLA o TRIANGULAR a partir de DXF
% =============================================================================
% Copyright (c) 2025 Nicolás Muñoz Guamán
% Email: nicolasjmunoz@gmail.com / nmunoz@mechaide.com
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

classdef mechaide_MESH
    methods (Static)
        % Generar malla de puntos a partir de un archivo DXF
        function mesh_data = generateMeshFromDXF(dxf_file, d, mesh_type, decimales, iterations, plot_contour, plot_base_mesh, plot_frontier, plot_mesh)
            fprintf('==================================================================\n');
            fprintf('               MECHAIDE Matlab v1.0 [21-08-2021]\n');
            fprintf('                    nmunoz@mechaide.com\n');
            fprintf('==================================================================\n');
            
            tic;
            plot_data = mechaide_DXF.plotSimple(dxf_file, d, false, false);
            fprintf('Importación de datos de DXF: %.4f s\n', toc);
            tic;
            puntos = mechaide_MESH.generatePointsByRegion(plot_data, d, decimales);
            data = mechaide_MESH.extractContoursAndHoles(puntos);
            fprintf('Generación de puntos en contornos: %.4f s\n', toc);
            tic;
            if plot_contour
                tic;
                mechaide_MESH.plotContour(data, true);
                fprintf('Gráfico de contorno: %.4f s\n', toc);
            end
            mesh_data = mechaide_MESH.generateRelaxedPointMesh(d, data.contorno, data.hoyos, iterations, mesh_type,  plot_base_mesh, plot_frontier);
            if plot_mesh
                tic;
                mechaide_MESH.plotFullMesh(mesh_data, false, true, true);
                fprintf('Gráfico de polygonos: %.4f s\n', toc);
            end
        end
        
        % Generar puntos a partir de datos DXF
        function newPoints = generatePointsByRegion(points_data, d, decimales)
            if nargin < 3
                decimales = 2;
            elseif nargin < 5
                decimales = 2;
            end
            
            % Obtener regiones únicas
            regionIds = unique(cellfun(@(p) p.region_id, points_data));
            newPoints = {};
            
            % Procesar cada región
            interpolate_edges = true;
            for r = 1:length(regionIds)
                regionId = regionIds(r);
                regionPoints = points_data(cellfun(@(p) p.region_id == regionId, points_data));
                
                regionNewPoints = mechaide_MESH.generateNewPoints(regionPoints, d, interpolate_edges);
                
                % Agregar puntos nuevos con id y region_id
                for i = 1:length(regionNewPoints)
                    newPoints{end+1} = struct(...
                        'id', i-1, ...
                        'region_id', regionId, ...
                        'x', mechaide_DXF.round(regionNewPoints{i}.x, decimales), ...
                        'y', mechaide_DXF.round(regionNewPoints{i}.y, decimales));
                end
            end
        end
        
        % Generar nuevos puntos a partir de un conjunto de puntos
        function newPoints = generateNewPoints(points, d, useExistingPoints)
            if nargin < 3
                useExistingPoints = true;
            end
            
            newPoints = {};
            for i = 1:length(points)
                currentPoint = points{i};
                nextPoint = points{mod(i, length(points)) + 1};
                
                % Calcular longitud del borde
                edgeLength = mechaide_MESH.distance(currentPoint, nextPoint);
                
                % Si la longitud es menor que d y useExistingPoints es true, conservar punto
                if edgeLength < d && useExistingPoints
                    newPoints{end+1} = currentPoint;
                end
                
                % Si la longitud es mayor o igual a d, generar nuevos puntos
                if edgeLength >= d
                    numNewPoints = floor(edgeLength / d);
                    for j = 0:numNewPoints
                        t = j / numNewPoints;
                        x = currentPoint.x + t * (nextPoint.x - currentPoint.x);
                        y = currentPoint.y + t * (nextPoint.y - currentPoint.y);
                        newPoints{end+1} = struct(...
                            'x', x, ...
                            'y', y, ...
                            'tipo', currentPoint.tipo, ...
                            'region_id', currentPoint.region_id, ...
                            'id', currentPoint.id, ...
                            'idu', currentPoint.idu);
                    end
                end
            end
        end
        
        % Calcular la distancia entre dos puntos
        function dist = distance(point1, point2)
            dist = sqrt((point2.x - point1.x)^2 + (point2.y - point1.y)^2);
        end
        
        % Extraer contorno y agujeros
        function data = extractContoursAndHoles(puntos)
            if isempty(puntos)
                data = struct('contorno', [], 'hoyos', {});
                return;
            end
            
            % Agrupar puntos por region_id
            regiones = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for i = 1:length(puntos)
                region_id = num2str(puntos{i}.region_id); % Convertir a string
                if ~regiones.isKey(region_id)
                    regiones(region_id) = zeros(0, 2);
                end
                regiones(region_id) = [regiones(region_id); [puntos{i}.x, puntos{i}.y]];
            end
            
            % Obtener region_id máximo como contorno
            keys = cellfun(@str2double, regiones.keys);
            if isempty(keys)
                data = struct('contorno', [], 'hoyos', {});
                return;
            end
            maxRegionId = num2str(max(keys)); % Convertir a string
            contorno = regiones(maxRegionId);
            
            % Resto de regiones como agujeros
            hoyos = {};
            for k = regiones.keys
                if ~strcmp(k{1}, maxRegionId)
                    hoyos{end+1} = regiones(k{1});
                end
            end
            
            data = struct('contorno', contorno, 'hoyos', {hoyos});
        end
        
        % Generar puntos relajados
        function data = generateRelaxedPointMesh(d, contorno, hoyos, iterations, mesh_type, plot_triangles, plot_frontier)
            if nargin < 4
                iterations = 1;
            end
            
            data = struct();
            data.centros = [];
            data.minX = [];
            data.minY = [];
            data.maxX = [];
            data.maxY = [];
            data.contorno = [];
            data.hoyos = {};
            data.iterations_points = {};
            
            tic;
            % Validar contorno
            if isempty(contorno) || size(contorno, 2) ~= 2 || size(contorno, 1) < 3
                warning('Contorno vacío o inválido. No se pueden generar puntos.');
                data.contorno = contorno;
                data.hoyos = hoyos;
                return;
            end
            
            % Verificar puntos únicos en el contorno
            unique_contorno = unique(contorno, 'rows');
            if size(unique_contorno, 1) < 3
                warning('El contorno es degenerado (menos de 3 puntos únicos). No se generan puntos.');
                data.contorno = contorno;
                data.hoyos = hoyos;
                return;
            end
            
            % No ordenar: asumir orden original del DXF para evitar desorden en concavidades
            if ~mechaide_MESH.isPolygonClosed(contorno)
                warning('El contorno no está cerrado. Cerrando manualmente.');
                contorno = [contorno; contorno(1, :)];
            end
            
            % Verificar área del contorno
            area = mechaide_MESH.polyarea(contorno);
            if area < 1e-6
                warning('El contorno es una línea recta (degenerado). No se pueden generar puntos dentro.');
                data.contorno = contorno;
                data.hoyos = hoyos;
                return;
            end
            
            % Validar agujeros
            valid_hoyos = {};
            for i = 1:length(hoyos)
                if ~isempty(hoyos{i}) && size(hoyos{i}, 2) == 2 && size(hoyos{i}, 1) >= 3
                    unique_hoyo = unique(hoyos{i}, 'rows');
                    if size(unique_hoyo, 1) >= 3
                        % No ordenar agujeros: asumir orden original
                        if ~mechaide_MESH.isPolygonClosed(hoyos{i})
                            hoyos{i} = [hoyos{i}; hoyos{i}(1, :)];
                        end
                        valid_hoyos{end+1} = hoyos{i};
                    else
                        warning('El hoyo %d es degenerado (menos de 3 puntos únicos). Ignorando.', i);
                    end
                end
            end
            hoyos = valid_hoyos;
            
            % Calcular límites
            limits = mechaide_MESH.findPointLimits(contorno);
            data.minX = limits(1);
            data.minY = limits(2);
            data.maxX = limits(3);
            data.maxY = limits(4);
            
            %fprintf('Límites: minX=%.2f, minY=%.2f, maxX=%.2f, maxY=%.2f\n', data.minX, data.minY, data.maxX, data.maxY);
            
            fprintf('Determinación de contorno y agujeros: %.4f s\n', toc);
            % Generar puntos con Poisson Disk Sampling
            tic;
            points = mechaide_MESH.generatePoissonDiskPoints(d, data.minX, data.minY, data.maxX, data.maxY, 30);
            %fprintf('Vértices generados por Poisson Disk: %d\n', size(points, 1));
            fprintf('Generación vertices por Poisson Disk: %.4f s\n', toc);
            
            % Filtrar puntos fuera del contorno o dentro de los agujeros
            tic;
            points = mechaide_MESH.excludePointsOutsideContourAndHoles(points, contorno, hoyos);
            %fprintf('Vértices después de excluir fuera de contorno/dentro de agujeros: %d\n', size(points, 1));
            
            if isempty(points)
                warning('No se encontraron puntos válidos dentro del contorno.');
                data.contorno = contorno;
                data.hoyos = hoyos;
                return;
            end
            
            % Definir puntos fijos (contorno + agujeros)
            fixedPoints = contorno;
            for i = 1:length(hoyos)
                if ~isempty(hoyos{i}) && size(hoyos{i}, 2) == 2
                    fixedPoints = [fixedPoints; hoyos{i}];
                end
            end
            fixedPoints = unique(fixedPoints, 'rows'); % Eliminar duplicados;
            
            % Filtrar puntos por distancia mínima a puntos fijos
            minDists = mechaide_MESH.minDistToPoints(points, fixedPoints);
            points = points(minDists >= d, :);
            %fprintf('Vértices después de filtrar distancia mínima a fijos: %d\n', size(points, 1));
            
            % Excluir puntos fijos de los puntos móviles
            points = points(~ismember(points, fixedPoints, 'rows'), :);
            fprintf('Filtrado de puntos no válidos: %.4f s\n', toc);
            fprintf('  - Vértices fijos (contorno + agujeros): %d\n', size(fixedPoints, 1))
            fprintf('  - Vértices móviles (interiores): %d\n', size(points, 1));
            
            fprintf('Iniciando teraciónes de relajación de Lloyd:\n');
            % Guardar puntos iniciales
            iterations_points = cell(1, iterations + 1);
            iterations_points{1} = points;
            % Aplicar relajación de Lloyd
            for i = 1:iterations
                tic; % Iniciar cronómetro
                points = mechaide_MESH.applyLloydRelaxationToMesh(points, fixedPoints, data.minX, data.minY, data.maxX, data.maxY, contorno, hoyos, d, mesh_type);
                fprintf('  #%d. Tiempo: %.4f s\n', i, toc);
                iterations_points{i+1} = points;
            end
            
            tic;
            % Combinar puntos móviles y fijos para la malla final
            centros = [points; fixedPoints];
            centros = unique(centros, 'rows'); % Eliminar duplicados
            data.centros = centros;
            data.coords = centros;
            
            % Generar conectividad según mesh_type
            switch lower(mesh_type)
                case 'triangle'
                    DT = delaunayTriangulation(centros);
                    connect = DT.ConnectivityList;
                    % Filtrar triángulos válidos (dentro de contorno, fuera de hoyos)
                    validTriangles = [];
                    for idx = 1:size(connect,1)
                        triangle = centros(connect(idx,:), :);
                        centroid = mean(triangle, 1);
                        inContorno = inpolygon(centroid(1), centroid(2), contorno(:,1), contorno(:,2));
                        inHoyo = false;
                        for j = 1:length(hoyos)
                            if inpolygon(centroid(1), centroid(2), hoyos{j}(:,1), hoyos{j}(:,2))
                                inHoyo = true;
                                break;
                            end
                        end
                        if inContorno && ~inHoyo
                            validTriangles = [validTriangles; connect(idx,:)];
                        end
                    end
                    data.connect = validTriangles;
                    data.he_data = mechaide_MESH.buildHalfEdgeMeshPolygon(centros, num2cell(validTriangles, 2));
                    
                    fprintf('Generación de malla triangular: %.4f s\n', toc);
                    fprintf('  - Triángulos validos: %d\n', size(validTriangles, 1));
                case 'dpolylla'
                    % Generar triangulación inicial
                    DT = delaunayTriangulation(centros);
                    connect = DT.ConnectivityList;
                    validTriangles = [];
                    for idx = 1:size(connect,1)
                        triangle = centros(connect(idx,:), :);
                        centroid = mean(triangle, 1);
                        inContorno = inpolygon(centroid(1), centroid(2), contorno(:,1), contorno(:,2));
                        inHoyo = false;
                        for j = 1:length(hoyos)
                            if inpolygon(centroid(1), centroid(2), hoyos{j}(:,1), hoyos{j}(:,2))
                                inHoyo = true;
                                break;
                            end
                        end
                        if inContorno && ~inHoyo
                            validTriangles = [validTriangles; connect(idx,:)];
                        end
                    end
                    
                    % Validar índices en validTriangles
                    if any(validTriangles(:) < 1 | validTriangles(:) > size(centros, 1))
                        error('Índices inválidos en validTriangles: min=%d, max=%d, num_vertices=%d', ...
                            min(validTriangles(:)), max(validTriangles(:)), size(centros, 1));
                    end
                    if any(isnan(centros(:)))
                        error('Valores NaN detectados en centros');
                    end
                    
                    fprintf('Generación de malla triangular: %.4f s\n', toc);
                    fprintf('  - Triángulos validos: %d\n', size(validTriangles, 1));
                    
                    tic;
                    % Gráfico 1: Malla original de triángulos (opcional)
                    if plot_triangles
                        figure; hold on;
                        triplot(validTriangles, centros(:,1), centros(:,2), 'k', 'LineWidth', 0.5);
                        plot(contorno(:,1), contorno(:,2), '-r', 'LineWidth', 1);
                        for j = 1:length(hoyos)
                            plot(hoyos{j}(:,1), hoyos{j}(:,2), '-m', 'LineWidth', 1);
                        end
                        title('Malla original de triángulos');
                        xlabel('X'); ylabel('Y'); axis equal;
                        hold off;
                    end
                    fprintf('Gráfico de malla de triángulos: %.4f s\n', toc);
                    
                    tic;
                    he_data = mechaide_MESH.buildHalfEdgeMeshPolygon(centros, num2cell(validTriangles, 2));
                    fprintf('Generación de estructura HE en triángulos: %.4f s\n', toc);
                    tic;
                    he_data = mechaide_MESH.labelMaxEdges(he_data);
                    fprintf('Identificación de aristas máximas: %.4f s\n', toc);
                    tic;
                    he_data = mechaide_MESH.labelFrontierEdges(he_data);
                    fprintf('Identificación de aristas de frontera: %.4f s\n', toc);
                    
                    tic;
                    frontierEdges = find([he_data.halfEdges.isFrontier]);
                    to_repair = mechaide_MESH.edgesToRepair(he_data);
                    fprintf('Identificación de aristas a reparar: %.4f s\n', toc);
                    fprintf('  - Aristas a reparar: %d aristas\n', length(to_repair));
                    
                    % Reparar y obtener nuevas aristas creadas
                    tic;
                    [frontierEdges, newEdges, he_data] = mechaide_MESH.repairBarrierEdgeTips(he_data, to_repair, frontierEdges);
                    fprintf('Reparación de aristas: %.4f s\n', toc);
                    
                    % Gráfico 2: Fronteras y reparaciones (opcional)
                    if plot_frontier
                        tic;
                        figure; hold on;
                        x_all = []; y_all = [];
                        for i = 1:length(frontierEdges)
                            he = he_data.halfEdges(frontierEdges(i));
                            if he.origin > 0 && he.origin <= length(he_data.vertices) && ...
                                    he.target > 0 && he.target <= length(he_data.vertices)
                                x_all = [x_all, he_data.vertices(he.origin).x, he_data.vertices(he.target).x, NaN];
                                y_all = [y_all, he_data.vertices(he.origin).y, he_data.vertices(he.target).y, NaN];
                            end
                        end
                        plot(x_all, y_all, '-b', 'LineWidth', 1);
                        
                        for i = 1:length(to_repair)
                            he = to_repair(i);
                            x = [he_data.vertices(he.origin).x, he_data.vertices(he.target).x];
                            y = [he_data.vertices(he.origin).y, he_data.vertices(he.target).y];
                            plot(x, y, '-r', 'LineWidth', 2);
                        end
                        
                        for i = 1:length(newEdges)
                            he = he_data.halfEdges(newEdges(i));
                            x = [he_data.vertices(he.origin).x, he_data.vertices(he.target).x];
                            y = [he_data.vertices(he.origin).y, he_data.vertices(he.target).y];
                            plot(x, y, '-g', 'LineWidth', 2);
                        end
                        
                        title('Aristas de frontera: azul = original, rojo = problemáticas, verde = reparadas');
                        xlabel('X'); ylabel('Y'); axis equal; hold off;
                        fprintf('Gráfico de aristas de frontera: %.4f s\n', toc);
                    end
                    
                    % Construir polígonos y conectividad
                    tic;
                    simple_polygons = mechaide_MESH.buildPolygons(he_data, frontierEdges);
                    connect_poly = mechaide_MESH.connectPolygons(simple_polygons, he_data);
                    meshData = mechaide_MESH.extractMeshDataPolygons(connect_poly);
                    he_data2 = mechaide_MESH.buildHalfEdgeMeshPolygon(meshData.coords, meshData.connect);
                    
                    % Resultado final
                    data.coords = meshData.coords;
                    data.connect = meshData.connect;
                    data.mesh = meshData.mesh;
                    data.he_data = he_data2;
                    
                    fprintf('Generación de malla polygonal: %.4f s\n', toc);
                    fprintf('  - Polígons validos: %d\n', size(simple_polygons, 1));
                    fprintf('  - Vértices totales: %d\n', size(data.centros, 1));
                otherwise
                    error('Tipo de malla no soportado: %s', mesh_type);
            end
            
            data.contorno = contorno;
            data.hoyos = hoyos;
            data.iterations_points = iterations_points;
            
            % Guardar en JSON
            tic;
            jsonStr = jsonencode(data, 'PrettyPrint', true);
            fid = fopen('malla.json', 'w');
            fprintf(fid, '%s', jsonStr);
            fclose(fid);
            fprintf('Exportando malla a archivo "malla.json": %.4f s\n', toc);
        end
        
        % Encuentra las aristas de frontera que deben ser reparadas
        function edges_to_repair = edgesToRepair(he_data)
            % 1. Filtrar todas las aristas que son fronteras
            allFrontierEdges = he_data.halfEdges([he_data.halfEdges.isFrontier]);
            
            % 2. Extraer todos los targets de las aristas de frontera
            targets = [allFrontierEdges.target];
            
            % 3. Contar ocurrencias de cada target
            [unique_targets, ~, idx] = unique(targets);
            targetCounts = accumarray(idx, 1);
            
            % 4. Filtrar aristas que NO son bordes y cuyo target aparece solo una vez
            edges_to_repair = allFrontierEdges(arrayfun(@(e) ...
                ~e.isBorder && targetCounts(unique_targets == e.target) == 1, ...
                allFrontierEdges));
        end
        
        % Construir polígonos y conectividad
        function polygons = buildPolygons(he_data, frontierEdges)
            nHE = length(he_data.halfEdges);
            
            % Preasignar vectores de acceso rápido
            HE_next   = zeros(1, nHE);
            HE_prev   = zeros(1, nHE);
            HE_origin = zeros(1, nHE);
            HE_target = zeros(1, nHE);
            HE_isFrontier = false(1, nHE);
            HE_twin   = zeros(1, nHE);
            
            for k = 1:nHE
                he = he_data.halfEdges(k);
                HE_next(k) = he.next;
                HE_prev(k) = he.prev;
                HE_origin(k) = he.origin;
                HE_target(k) = he.target;
                HE_isFrontier(k) = he.isFrontier;
                if isempty(he.twin)
                    HE_twin(k) = 0;
                else
                    HE_twin(k) = he.twin;
                end
            end
            
            visitedEdges = false(1, nHE);
            polygons = {};
            
            %fprintf('Iniciando construcción de polígonos...');
            
            for k = 1:length(frontierEdges)
                edgeIdx = frontierEdges(k);
                if visitedEdges(edgeIdx)
                    continue;
                end
                
                %fprintf('➤ Nuevo polígono iniciado desde arista %d\n', edgeIdx);
                
                polygon = [];
                vertices = [];
                actual = edgeIdx;
                
                % antes de entrar al while, inicializa localVisited
                localVisited = false(1, nHE);  % nHE ya definido fuera
                
                while true
                    % protección: si actual ya fue visto en este polígono -> romper (evita bucle infinito)
                    if localVisited(actual)
                        %warning('buildPolygons: ciclo detectado dentro del mismo polígono empezando en arista %d. Rompiendo.', edgeIdx);
                        break;
                    end
                    localVisited(actual) = true;
                    
                    % Triángulo actual
                    triEdges = [actual, HE_next(actual), HE_prev(actual)];
                    
                    % Separar aristas de frontera y no-frontera
                    isFront = HE_isFrontier(triEdges);
                    seedTriangleEdges = triEdges(isFront);
                    emptyEdges = triEdges(~isFront);
                    
                    % Agregar aristas de frontera al polígono
                    for e = seedTriangleEdges
                        if ~visitedEdges(e)
                            visitedEdges(e) = true;
                            polygon(end+1) = e;
                            vertices(end+1:end+2) = [HE_origin(e), HE_target(e)];
                            %fprintf('   + Arista %d agregada (origen %d → target %d)\n', ...
                            %    e, HE_origin(e), HE_target(e));
                        end
                    end
                    
                    % Mostrar estado actual del polígono
                    %fprintf('   Polígono actual: [');
                    %fprintf('%d ', polygon);
                    %fprintf(']\n');
                    
                    % Verificar cierre del polígono
                    if numel(seedTriangleEdges) == 3
                        %fprintf('   ✓ Polígono cerrado (3 aristas de frontera)\n');
                        break;
                    end
                    
                    vertCounts = histc(vertices, 1:length(he_data.vertices));
                    if all(vertCounts(vertices) == 2)
                        %fprintf('   ✓ Polígono cerrado por vértices\n');
                        break;
                    end
                    
                    % Avanzar al twin de la última arista no-frontera
                    if ~isempty(emptyEdges)
                        twinIdx = HE_twin(emptyEdges(end));
                        if twinIdx > 0 && ~HE_isFrontier(twinIdx) && ~visitedEdges(twinIdx)
                            actual = twinIdx;
                            continue;
                        end
                    end
                    break;
                end
                
                %fprintf('✔ Polígono finalizado con %d aristas\n\n', length(polygon));
                polygons{end+1} = polygon;
            end
            
        end
        
        % Reparar aristas con vértices problemáticos
        function [frontierEdges, newEdges, he_data] = repairBarrierEdgeTips(he_data, to_repair, frontierEdges)
            newEdges = [];
            
            for ii = 1:length(to_repair)
                problematicEdge = to_repair(ii);
                
                % --- Obtener índice del vértice problemático (asegurarlo como int) ---
                if isstruct(problematicEdge.target)
                    problematicVertex = problematicEdge.target.index;
                else
                    problematicVertex = problematicEdge.target;
                end
                problematicVertex = double(problematicVertex);
                
                % --- Buscar aristas NO frontera conectadas al vértice problemático ---
                candidates = [];
                for j = 1:length(he_data.halfEdges)
                    he = he_data.halfEdges(j);
                    if ~he.isFrontier
                        he_origin = he.origin; if isstruct(he_origin), he_origin = he_origin.index; end
                        he_target = he.target; if isstruct(he_target), he_target = he_target.index; end
                        he_origin = double(he_origin); he_target = double(he_target);
                        
                        if he_origin == problematicVertex || he_target == problematicVertex
                            candidates(end+1) = j; %#ok<AGROW>
                        end
                    end
                end
                
                %fprintf('Caso %d: vertex %d, candidatos encontrados: %d\n', ii, problematicVertex, length(candidates));
                if isempty(candidates)
                    %warning('  -> No hay candidatos no-frontier conectados al vértice %d. Saltando.', problematicVertex);
                    continue;
                end
                
                % --- Escoger el candidato más largo (edgeLengthSquared devuelve número) ---
                longestIdx = candidates(1);
                maxL = mechaide_MESH.edgeLengthSquared(he_data.halfEdges(longestIdx), he_data);
                for k = 2:length(candidates)
                    L = mechaide_MESH.edgeLengthSquared(he_data.halfEdges(candidates(k)), he_data);
                    if L > maxL
                        maxL = L;
                        longestIdx = candidates(k);
                    end
                end
                
                % --- Log de candidato seleccionado ---
                heL = he_data.halfEdges(longestIdx);
                % normalizar origin/target a índices numéricos
                oL = heL.origin; if isstruct(oL), oL = oL.index; end
                tL = heL.target; if isstruct(tL), tL = tL.index; end
                %fprintf('  -> Seleccionada arista %d (or %d -> tgt %d), longitud^2=%.4g\n', longestIdx, oL, tL, maxL);
                
                % --- Promocionar la arista y su twin, SI AÚN NO ESTÁN en frontierEdges ---
                if ~ismember(longestIdx, frontierEdges)
                    frontierEdges(end+1) = longestIdx;            % append
                    newEdges(end+1) = longestIdx;                 % append
                    he_data.halfEdges(longestIdx).isFrontier = true;
                    %fprintf('     + Promovida %d\n', longestIdx);
                else
                    %fprintf('     - %d ya era frontier\n', longestIdx);
                end
                
                twin = he_data.halfEdges(longestIdx).twin;
                if ~isempty(twin) && isnumeric(twin)
                    if ~ismember(twin, frontierEdges)
                        frontierEdges(end+1) = twin;
                        newEdges(end+1) = twin;
                        he_data.halfEdges(twin).isFrontier = true;
                        %fprintf('     + Promovida twin %d\n', twin);
                    else
                        %fprintf('     - twin %d ya era frontier\n', twin);
                    end
                else
                    %fprintf('     - Twin inexistente o inválido para %d\n', longestIdx);
                end
                
            end % for to_repair
            
            % Asegurar unicidad en frontierEdges y newEdges (mantener orden original)
            frontierEdges = unique(frontierEdges, 'stable');
            newEdges = unique(newEdges, 'stable');
            
            %fprintf('Reparación finalizada. Nuevas aristas: %d', length(newEdges));
        end
        
        
        function length_sq = edgeLengthSquared(he, he_data)
            % Obtener structs de vértices usando índices
            origin = he_data.vertices(he.origin); % Acceder al struct del vértice origen
            target = he_data.vertices(he_data.halfEdges(he.next).origin); % Acceder al struct del vértice destino
            dx = target.x - origin.x;
            dy = target.y - origin.y;
            length_sq = dx^2 + dy^2;
        end
        
        
        function he_data = labelMaxEdges(he_data)
            for i = 1:length(he_data.faces)
                face = he_data.faces(i);
                heIndices = [face.edge];
                heIndices(2) = he_data.halfEdges(heIndices(1)).next;
                heIndices(3) = he_data.halfEdges(heIndices(2)).next;
                lengths = zeros(1, 3);
                for j = 1:3
                    he = he_data.halfEdges(heIndices(j));
                    lengths(j) = mechaide_MESH.edgeLengthSquared(he, he_data); % Usar la función corregida
                end
                [~, maxIdx] = max(lengths);
                he_data.halfEdges(heIndices(maxIdx)).isMax = true;
            end
        end
        
        function he_data = labelFrontierEdges(he_data)
            frontierEdges = [];
            for i = 1:length(he_data.halfEdges)
                he = he_data.halfEdges(i);
                he_data.halfEdges(i).isFrontier = false;
                if he.isBorder
                    he_data.halfEdges(i).isFrontier = true;
                    frontierEdges(end+1) = i;
                elseif ~isempty(he.twin) && he.twin > 0 && he.twin <= length(he_data.halfEdges)
                    twinHe = he_data.halfEdges(he.twin);
                    if ~he.isMax && ~twinHe.isMax
                        he_data.halfEdges(i).isFrontier = true;
                        frontierEdges(end+1) = i;
                    end
                    if he.origin == twinHe.target && he.target == twinHe.origin
                        continue;
                    else
                        warning('Twin inconsistente en he %d: origin=%d, target=%d; twin origin=%d, target=%d', ...
                            i, he.origin, he.target, twinHe.origin, twinHe.target);
                    end
                end
            end
            fprintf('labelFrontierEdges: %d aristas de frontera\n', length(frontierEdges));
        end
        
        function he_data = buildHalfEdgeMeshPolygon(coords, connect)
            % Inicializar estructuras
            vertices = struct('x', num2cell(coords(:,1)), 'y', num2cell(coords(:,2)), ...
                'index', num2cell(1:size(coords,1))', 'incidentEdge', [], 'isBorder', false);
            halfEdges = struct('origin', {}, 'target', {}, 'twin', {}, 'next', {}, 'prev', {}, ...
                'face', {}, 'key', {}, 'isBorder', {}, 'isMax', {}, 'isFrontier', {}, 'marked', {}, 'index', {});
            faces = struct('edge', {}, 'index', {});
            edgeMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
            heIndex = 1;
            faceIndex = 1;
            
            % Procesar cada polígono
            for i = 1:length(connect)
                polygon = connect{i};
                if length(polygon) < 3 || any(polygon < 1 | polygon > size(coords,1))
                    continue;
                end
                face = struct('edge', [], 'index', faceIndex);
                faces(faceIndex) = face;
                faceHalfEdges = [];
                
                % Crear half-edges para cada lado del polígono
                for j = 1:length(polygon)
                    originIdx = polygon(j);
                    targetIdx = polygon(mod(j, length(polygon)) + 1);
                    if originIdx == targetIdx
                        %warning('Arista degenerada en polígono %d: origen=%d, destino=%d', i, originIdx, targetIdx);
                        continue;
                    end
                    % Crear half-edge
                    he = struct('origin', originIdx, ... % Índice en lugar de struct
                        'target', targetIdx, ... % Índice en lugar de struct
                        'twin', [], ...
                        'next', [], ...
                        'prev', [], ...
                        'face', faceIndex, ...
                        'key', sprintf('%d-%d', originIdx, targetIdx), ...
                        'isBorder', false, ...
                        'isMax', false, ...
                        'isFrontier', false, ...
                        'marked', false, ...
                        'index', heIndex);
                    halfEdges(heIndex) = he;
                    faceHalfEdges(end+1) = heIndex;
                    
                    % Asignar incidentEdge
                    if isempty(vertices(originIdx).incidentEdge)
                        vertices(originIdx).incidentEdge = heIndex;
                    end
                    
                    % Manejar twins
                    key = sprintf('%d-%d', min(originIdx, targetIdx), max(originIdx, targetIdx));
                    if edgeMap.isKey(key)
                        twinIdx = edgeMap(key);
                        halfEdges(heIndex).twin = twinIdx;
                        halfEdges(twinIdx).twin = heIndex;
                        edgeMap.remove(key);
                    else
                        edgeMap(key) = heIndex;
                    end
                    heIndex = heIndex + 1;
                end
                
                % Asignar next y prev
                for j = 1:length(faceHalfEdges)
                    halfEdges(faceHalfEdges(j)).next = faceHalfEdges(mod(j, length(faceHalfEdges)) + 1);
                    halfEdges(faceHalfEdges(j)).prev = faceHalfEdges(mod(j-2, length(faceHalfEdges)) + 1);
                end
                faces(faceIndex).edge = faceHalfEdges(1);
                faceIndex = faceIndex + 1;
            end
            
            % Marcar aristas de borde
            for key = edgeMap.keys
                heIdx = edgeMap(key{1});
                halfEdges(heIdx).isBorder = true;
                halfEdges(heIdx).twin = [];
                vertices(halfEdges(heIdx).origin).isBorder = true;
            end
            
            % Limpiar estructuras
            halfEdges = halfEdges(1:heIndex-1);
            faces = faces(1:faceIndex-1);
            
            % Retornar he_data
            he_data = struct('vertices', vertices, 'halfEdges', halfEdges, 'faces', faces);
        end
        
        function meshData = extractMeshDataPolygons(polygons)
            coordsMap = containers.Map('KeyType','char','ValueType','double');
            coords = [];
            connect = {};
            idx = 1;
            
            for i = 1:length(polygons)
                polygon = polygons{i};
                element = [];
                previousIndex = -1;
                for j = 1:length(polygon)
                    point = polygon{j};
                    key = sprintf('%.12f,%.12f', point(1), point(2));
                    if ~coordsMap.isKey(key)
                        coordsMap(key) = idx;
                        coords(idx,:) = point;
                        idx = idx + 1;
                    end
                    pointIndex = coordsMap(key);
                    if pointIndex ~= previousIndex
                        element(end+1) = pointIndex;
                    end
                    previousIndex = pointIndex;
                end
                if length(element) >= 3
                    connect{end+1} = element;
                else
                    warning('Polígono %d con menos de 3 puntos, omitido', i);
                end
            end
            
            meshData.coords = coords;
            meshData.connect = connect;
            meshData.mesh = polygons;
        end
        function connectedPolygons = connectPolygons(polygons, he_data)
            % polygons: cell array donde cada elemento es vector de half-edge indices
            % he_data: estructura original con he_data.halfEdges y he_data.vertices
            connectedPolygons = cell(1, length(polygons));
            for i = 1:length(polygons)
                connectedPolygons{i} = mechaide_MESH.ordenarEdges(polygons{i}, he_data);
            end
        end
        function vertexCoordsCell = ordenarEdges(edges, he_data)
            % edges: puede ser
            %   - vector de índices (ej. [1351 1352 ...])  <-- lo que buildPolygons produce, o
            %   - array de structs con campos origin/target (menos probable aquí)
            %
            % he_data: necesario para acceder a he_data.halfEdges y he_data.vertices
            %
            % Devuelve: cell array 1 x N de puntos [x,y] para ese polígono:
            %   vertexCoordsCell = { [x1,y1], [x2,y2], ... }
            
            vertexCoordsCell = {};
            if isempty(edges)
                return;
            end
            
            % Si edges son structs, convertir a índices (detectar)
            if isstruct(edges)
                % convertir a índices si el struct tiene campo 'index'
                if isfield(edges, 'index')
                    edgeIdxs = [edges.index];
                else
                    error('ordenarEdges: structs sin campo index.');
                end
            else
                edgeIdxs = edges; % asumo vector numérico de índices
            end
            
            % Obtener array de half-edge structs
            HE = he_data.halfEdges;
            
            % Empezar con la primer arista del conjunto
            firstEdgeIdx = edgeIdxs(1);
            if firstEdgeIdx < 1 || firstEdgeIdx > length(HE)
                warning('ordenarEdges: primer índice fuera de rango: %d', firstEdgeIdx);
                return;
            end
            edgeActualIdx = firstEdgeIdx;
            
            % primer vértice (origin)
            originIdx = HE(edgeActualIdx).origin;
            % normalizar si origin es struct
            if isstruct(originIdx) && isfield(originIdx, 'index')
                originIdx = originIdx.index;
            end
            primerIndex = originIdx;
            
            % agregar primer punto
            v = he_data.vertices(primerIndex);
            vertexCoordsCell{end+1} = [v.x, v.y];
            
            maxSteps = numel(edgeIdxs) * 4 + 50; % protección
            steps = 0;
            
            while true
                steps = steps + 1;
                if steps > maxSteps
                    warning('ordenarEdges: posible bucle infinito, rompiendo (edge %d).', firstEdgeIdx);
                    break;
                end
                
                % añadir punto target del edgeActual
                targ = HE(edgeActualIdx).target;
                if isstruct(targ) && isfield(targ, 'index'), targ = targ.index; end
                v = he_data.vertices(targ);
                vertexCoordsCell{end+1} = [v.x, v.y];
                
                % buscar siguiente edge dentro del conjunto edgeIdxs cuyo origin == targ
                % construir lista de origins de edgeIdxs una sola vez:
                % (notar: small set, find es aceptable)
                siguienteEdgePos = find(arrayfun(@(e) ...
                    get_origin_index(HE(e)), edgeIdxs) == targ, 1);
                
                if isempty(siguienteEdgePos)
                    break;
                end
                
                siguienteEdgeIdx = edgeIdxs(siguienteEdgePos);
                % si el siguiente edge cierra ciclo con primerIndex -> salir
                nextTarget = HE(siguienteEdgeIdx).target;
                if isstruct(nextTarget) && isfield(nextTarget, 'index'), nextTarget = nextTarget.index; end
                if nextTarget == primerIndex
                    % añadir el target final (que es primerIndex) y terminar
                    v = he_data.vertices(nextTarget);
                    vertexCoordsCell{end+1} = [v.x, v.y];
                    break;
                end
                
                % avanzar
                edgeActualIdx = siguienteEdgeIdx;
            end
            % helper function (local)
            function oidx = get_origin_index(he)
                o = he.origin;
                if isstruct(o) && isfield(o,'index')
                    oidx = o.index;
                else
                    oidx = o;
                end
            end
        end
        function area = polyarea(polygon)
            n = size(polygon, 1);
            area = 0;
            for i = 1:n
                x1 = polygon(i, 1);
                y1 = polygon(i, 2);
                x2 = polygon(mod(i, n) + 1, 1);
                y2 = polygon(mod(i, n) + 1, 2);
                area = area + (x1 * y2 - x2 * y1);
            end
            area = abs(area) / 2;
        end
        function minDists = minDistToPoints(points, fixedPoints)
            n = size(points, 1);
            minDists = inf(n, 1);
            for i = 1:n
                for j = 1:size(fixedPoints, 1)
                    dx = points(i, 1) - fixedPoints(j, 1);
                    dy = points(i, 2) - fixedPoints(j, 2);
                    dist = sqrt(dx^2 + dy^2);
                    if dist < minDists(i)
                        minDists(i) = dist;
                    end
                end
            end
        end
        function limits = findPointLimits(points)
            % Calcula los límites de un conjunto de puntos
            % Entrada:
            %   points: Matriz [x, y] de coordenadas
            % Salida:
            %   Vector [minX, minY, maxX, maxY]
            minX = min(points(:, 1));
            minY = min(points(:, 2));
            maxX = max(points(:, 1));
            maxY = max(points(:, 2));
            limits = [minX, minY, maxX, maxY];
        end
        function points = generatePoissonDiskPoints(d, minX, minY, maxX, maxY, iterations)
            % Genera puntos usando Poisson Disk Sampling en un rectángulo definido
            % d: distancia mínima entre puntos
            % minX, minY, maxX, maxY: límites del dominio
            % iterations: intentos por punto
            
            % Tamaño del dominio
            width = maxX - minX;
            height = maxY - minY;
            
            % Tamaño de celda para la grilla auxiliar
            cellSize = d / sqrt(2);
            
            % Grilla de referencia (índices de puntos)
            gridWidth = ceil(width / cellSize);
            gridHeight = ceil(height / cellSize);
            grid = -ones(gridHeight, gridWidth);
            
            % Lista de puntos y de activos
            points = zeros(0, 2);
            active = zeros(0, 2);
            
            % Punto inicial aleatorio
            initialPoint = [minX + rand() * width, minY + rand() * height];
            points(1, :) = initialPoint;
            active(1, :) = initialPoint;
            
            % Guardar en la grilla
            gx = floor((initialPoint(1) - minX) / cellSize) + 1;
            gy = floor((initialPoint(2) - minY) / cellSize) + 1;
            grid(gy, gx) = 1;
            
            while ~isempty(active)
                idx = randi(size(active, 1));  % punto activo aleatorio
                point = active(idx, :);
                found = false;
                
                for i = 1:iterations
                    % Generar nuevo punto a una distancia entre d y 2d
                    r = d * (1 + rand());
                    theta = 2 * pi * rand();
                    newPoint = point + r * [cos(theta), sin(theta)];
                    
                    % Verificar que esté dentro de límites
                    if newPoint(1) >= minX && newPoint(1) <= maxX && ...
                            newPoint(2) >= minY && newPoint(2) <= maxY
                        
                        % Verificar distancia mínima usando la grilla
                        gx = floor((newPoint(1) - minX) / cellSize) + 1;
                        gy = floor((newPoint(2) - minY) / cellSize) + 1;
                        
                        % Revisar vecinos
                        ok = true;
                        for ix = max(gx-2, 1):min(gx+2, gridWidth)
                            for iy = max(gy-2, 1):min(gy+2, gridHeight)
                                pIndex = grid(iy, ix);
                                if pIndex ~= -1
                                    dist2 = sum((points(pIndex,:) - newPoint).^2);
                                    if dist2 < d^2
                                        ok = false;
                                        break;
                                    end
                                end
                            end
                            if ~ok, break; end
                        end
                        
                        % Si es válido, añadirlo
                        if ok
                            points(end+1, :) = newPoint;
                            active(end+1, :) = newPoint;
                            grid(gy, gx) = size(points, 1);
                            found = true;
                            break;
                        end
                    end
                end
                
                % Si no se encontró, eliminar de activos
                if ~found
                    active(idx, :) = [];
                end
            end
        end
        function points = excludePointsOutsideContourAndHoles(points, contorno, hoyos)
            if isempty(points) || isempty(contorno) || size(contorno, 1) < 3
                %warning('Vértices o contorno vacíos/inválidos. Retornando vacío.');
                points = [];
                return;
            end
            
            % Verificar si el contorno es válido
            unique_contorno = unique(contorno, 'rows');
            if size(unique_contorno, 1) < 3
                warning('El contorno es degenerado (menos de 3 puntos únicos). Retornando vacío.');
                points = [];
                return;
            end
            
            % Verificar autointersecciones
            if mechaide_MESH.hasSelfIntersections(contorno)
                %warning('El contorno tiene autointersecciones, lo que puede causar errores en inpolygon.');
            end
            
            % Asegurar que el contorno esté en sentido antihorario
            area = mechaide_MESH.polyarea(contorno);
            if area < 0
                contorno = flipud(contorno);
                fprintf('Contorno invertido a sentido antihorario.\n');
            end
            
            % Verificar si el contorno está cerrado
            if ~mechaide_MESH.isPolygonClosed(contorno)
                warning('El contorno no está cerrado. Cerrando manualmente.');
                contorno = [contorno; contorno(1, :)];
            end
            
            % Verificar área
            area = mechaide_MESH.polyarea(contorno);
            if area < 1e-6
                warning('El contorno tiene área cercana a cero (degenerado). No se generan puntos.');
                points = [];
                return;
            end
            
            % Eliminar puntos duplicados consecutivos en el contorno
            contorno = unique(contorno, 'rows', 'stable');
            if ~isequal(contorno(1, :), contorno(end, :))
                contorno = [contorno; contorno(1, :)];
            end
            
            % Filtrar puntos dentro del contorno
            [in, on] = inpolygon(points(:,1), points(:,2), contorno(:,1), contorno(:,2));
            valid_points = in & ~on;
            if ~any(valid_points)
                warning('Ningún punto está estrictamente dentro del contorno. Verifica la geometría.');
                figure;
                plot(contorno(:,1), contorno(:,2), 'b-', 'LineWidth', 1);
                hold on;
                plot(points(:,1), points(:,2), 'r.', 'MarkerSize', 10);
                plot(points(valid_points,1), points(valid_points,2), 'g.', 'MarkerSize', 10);
                axis equal;
                title('Depuración: Puntos (rojo), Contorno (azul), Válidos (verde)');
                hold off;
            end
            discarded_points = points(~valid_points, :);
            if ~isempty(discarded_points)
                %fprintf('Primeros 5 puntos descartados:\n');
                %disp(discarded_points(1:min(5, size(discarded_points, 1)), :));
            end
            points = points(valid_points, :);
            
            % Filtrar puntos dentro de los agujeros
            for i = 1:length(hoyos)
                if ~isempty(hoyos{i}) && size(hoyos{i}, 2) == 2 && size(hoyos{i}, 1) >= 3
                    unique_hoyo = unique(hoyos{i}, 'rows');
                    if size(unique_hoyo, 1) < 3
                        warning('El hoyo %d es degenerado (menos de 3 puntos únicos). Ignorando.', i);
                        continue;
                    end
                    if ~mechaide_MESH.isPolygonClosed(hoyos{i})
                        hoyos{i} = [hoyos{i}; hoyos{i}(1, :)];
                    end
                    area_hoyo = mechaide_MESH.polyarea(hoyos{i});
                    if area_hoyo > 0
                        hoyos{i} = flipud(hoyos{i});
                    end
                    if mechaide_MESH.hasSelfIntersections(hoyos{i})
                    end
                    [in, on] = inpolygon(points(:,1), points(:,2), hoyos{i}(:,1), hoyos{i}(:,2));
                    points = points(~(in | on), :);
                end
            end
            
            % Filtrado adicional para puntos cercanos al contorno
            if ~isempty(points)
                min_dist = zeros(size(points, 1), 1);
                for j = 1:size(points, 1)
                    min_dist(j) = mechaide_MESH.distanceToPolygon(points(j, :), contorno);
                end
                valid_dist = min_dist > 1e-3; % Ajustar umbral según la escala
                points = points(valid_dist, :);
            end
            
        end
        function intersect = segmentsIntersect(p1, p2, q1, q2)
            % Verifica si dos segmentos de línea se intersectan
            % Basado en el algoritmo de orientación y bounding box
            o1 = mechaide_MESH.orientation(p1, p2, q1);
            o2 = mechaide_MESH.orientation(p1, p2, q2);
            o3 = mechaide_MESH.orientation(q1, q2, p1);
            o4 = mechaide_MESH.orientation(q1, q2, p2);
            
            % Caso general: diferentes orientaciones
            intersect = (o1 ~= o2 && o3 ~= o4);
            
            % Casos especiales: colineales
            if o1 == 0 && o2 == 0 && o3 == 0 && o4 == 0
                % Verificar si los segmentos se superponen
                intersect = mechaide_MESH.onSegment(p1, q1, p2) || mechaide_MESH.onSegment(p1, q2, p2) || ...
                    mechaide_MESH.onSegment(q1, p1, q2) || mechaide_MESH.onSegment(q1, p2, q2);
            end
        end
        function o = orientation(p, q, r)
            % Calcula la orientación de tres puntos (p, q, r)
            % 0: colineales, 1: horario, -1: antihorario
            val = (q(2) - p(2)) * (r(1) - q(1)) - (q(1) - p(1)) * (r(2) - q(2));
            if abs(val) < 1e-10
                o = 0; % Colineal
            elseif val > 0
                o = 1; % Horario
            else
                o = -1; % Antihorario
            end
        end
        function on_seg = onSegment(p, q, r)
            % Verifica si el punto q está en el segmento pr
            on_seg = q(1) <= max(p(1), r(1)) && q(1) >= min(p(1), r(1)) && ...
                q(2) <= max(p(2), r(2)) && q(2) >= min(p(2), r(2));
        end
        function has_intersections = hasSelfIntersections(polygon)
            % Verifica si un polígono tiene autointersecciones
            has_intersections = false;
            n = size(polygon, 1);
            if n < 4
                return;
            end
            
            % Comparar cada segmento con todos los demás
            for i = 1:n-1
                p1 = polygon(i, :);
                p2 = polygon(mod(i, n) + 1, :);
                for j = i+2:n-1
                    % Evitar segmentos adyacentes
                    if j == mod(i, n) + 1 || j == i-1
                        continue;
                    end
                    q1 = polygon(j, :);
                    q2 = polygon(mod(j, n) + 1, :);
                    
                    % Verificar intersección entre segmentos [p1, p2] y [q1, q2]
                    if mechaide_MESH.segmentsIntersect(p1, p2, q1, q2)
                        has_intersections = true;
                        return;
                    end
                end
            end
        end
        function closed = isPolygonClosed(polygon)
            % Verifica si un polígono está cerrado
            % Entrada:
            %   polygon: Matriz [x, y] de puntos
            % Salida:
            %   closed: True si el primer y último punto son iguales
            closed = isequal(polygon(1, :), polygon(end, :));
        end
        function dist = distanceToPolygon(point, polygon)
            % Calcula la distancia mínima de un punto a un polígono
            % Entrada:
            %   point: Vector [x, y]
            %   polygon: Matriz [x, y] de puntos del polígono
            % Salida:
            %   dist: Distancia mínima
            dist = inf;
            for i = 1:size(polygon, 1)
                a = polygon(i, :);
                b = polygon(mod(i, size(polygon, 1)) + 1, :);
                dist = min(dist, mechaide_MESH.pointToSegmentDistance(point, a, b));
            end
        end
        function dist = pointToSegmentDistance(p, a, b)
            % Calcula la distancia de un punto a un segmento de línea
            % Entrada:
            %   p: Punto [x, y]
            %   a, b: Puntos del segmento [x, y]
            % Salida:
            %   dist: Distancia al segmento
            px = p(1); py = p(2);
            ax = a(1); ay = a(2);
            bx = b(1); by = b(2);
            
            dx = bx - ax;
            dy = by - ay;
            
            if dx == 0 && dy == 0
                % a y b son el mismo punto
                dist = sqrt((px - ax)^2 + (py - ay)^2);
                return;
            end
            
            % Calcular parámetro t
            t = ((px - ax) * dx + (py - ay) * dy) / (dx^2 + dy^2);
            
            if t < 0
                closest = a;
            elseif t > 1
                closest = b;
            else
                closest = [ax + t * dx, ay + t * dy];
            end
            
            dist = sqrt((px - closest(1))^2 + (py - closest(2))^2);
        end
        function relaxedPoints = applyLloydRelaxationToMesh(points, fixedPoints, minX, minY, maxX, maxY, contorno, hoyos, d, mesh_type)
            % Aplica relajación de Lloyd a puntos móviles
            % Entrada:
            %   points: Matriz [x, y] de puntos móviles
            %   fixedPoints: Matriz [x, y] de puntos fijos
            %   minX, minY, maxX, maxY: Límites del área
            %   contorno: Matriz [x, y] del contorno
            %   hoyos: Cell array de matrices [x, y] de agujeros
            %   d: Distancia mínima
            %   mesh_type: 'triangle' o 'polygon'
            % Salida:
            %   relaxedPoints: Matriz [x, y] de puntos relajados
            
            % Combinar puntos móviles y fijos, eliminando duplicados
            combinedPoints = unique([points; fixedPoints], 'rows', 'stable');
            
            % Calcular triangulación de Delaunay y diagrama de Voronoi
            DT = delaunayTriangulation(combinedPoints);
            [V, R] = voronoiDiagram(DT);
            
            % Crear conjunto de puntos fijos para búsqueda rápida
            fixedPointsSet = fixedPoints;
            
            % Asegurar que contorno y agujeros estén cerrados
            if ~mechaide_MESH.isPolygonClosed(contorno)
                contorno = [contorno; contorno(1, :)];
            end
            for i = 1:length(hoyos)
                if ~mechaide_MESH.isPolygonClosed(hoyos{i})
                    hoyos{i} = [hoyos{i}; hoyos{i}(1, :)];
                end
            end
            
            % Aplicar relajación de Lloyd a puntos móviles
            relaxedPoints = points;
            for i = 1:size(points, 1)
                point = points(i, :);
                if ismember(point, fixedPointsSet, 'rows')
                    continue; % Saltar puntos fijos
                end
                
                % Obtener celda de Voronoi
                idx = find(all(abs(combinedPoints - point) < 1e-6, 2), 1); % Buscar índice del punto
                if isempty(idx) || isempty(R{idx}) || any(any(isinf(V(R{idx}, :))))
                    continue; % Celdas inválidas
                end
                
                cellVertices = V(R{idx}, :);
                if isempty(cellVertices) || any(any(isinf(cellVertices)))
                    continue; % Celdas inválidas
                end
                
                % Calcular centroide de la celda
                centroid = mean(cellVertices, 1);
                if any(isnan(centroid))
                    continue;
                end
                
                % Verificar si el centroide está dentro del contorno
                %if ~inpolygon(centroid(1), centroid(2), contorno(:, 1), contorno(:, 2))
                %    continue;
                %end
                
                % Verificar si el centroide está dentro de algún agujero
                %in_hoyo = false;
                %for j = 1:length(hoyos)
                %    if inpolygon(centroid(1), centroid(2), hoyos{j}(:, 1), hoyos{j}(:, 2))
                %        in_hoyo = true;
                %        break;
                %    end
                %end
                %if in_hoyo
                %    continue;
                %end
                
                % Verificar distancia mínima al contorno
                distContorno = mechaide_MESH.distanceToPolygon(centroid, contorno);
                if strcmp(mesh_type, 'polygon')
                    dist_threshold = d / 3;
                else
                    dist_threshold = d;
                end
                if distContorno < dist_threshold
                    continue;
                end
                
                % Verificar distancia mínima a los agujeros
                for j = 1:length(hoyos)
                    distHoyo = mechaide_MESH.distanceToPolygon(centroid, hoyos{j});
                    if strcmp(mesh_type, 'polygon')
                        dist_threshold = d * 3;
                    else
                        dist_threshold = d;
                    end
                    if distHoyo < dist_threshold
                        continue;
                    end
                end
                
                % Verificar distancia mínima a puntos fijos
                distPuntoFijo = min(sqrt(sum((fixedPoints - centroid).^2, 2)));
                if distPuntoFijo < d
                    continue;
                end
                
                % Verificar distancia mínima a otros puntos móviles
                otherPoints = points(setdiff(1:size(points, 1), i), :);
                if ~isempty(otherPoints)
                    distOtroPunto = min(sqrt(sum((otherPoints - centroid).^2, 2)));
                    if distOtroPunto < d
                        continue;
                    end
                end
                
                % Mover punto al centroide
                relaxedPoints(i, :) = centroid;
            end
        end
        function plotContour(data, enable)
            if nargin < 2
                enable = true;
            end
            if ~enable
                return;
            end
            
            if ~isstruct(data) || ~isfield(data, 'contorno')
                warning('Estructura data inválida o sin contorno.');
                return;
            end
            
            figure;
            hold on;
            if ~isempty(data.contorno) && size(data.contorno, 2) == 2
                plot(data.contorno(:,1), data.contorno(:,2), 'b-', 'LineWidth', 1);
            else
                warning('Contorno inválido o vacío.');
            end
            
            if isfield(data, 'hoyos') && ~isempty(data.hoyos)
                for i = 1:length(data.hoyos)
                    if size(data.hoyos{i}, 2) == 2
                        plot(data.hoyos{i}(:,1), data.hoyos{i}(:,2), 'r-', 'LineWidth', 1);
                    end
                end
            end
            hold off;
            axis equal;
            title('Contorno');
        end
        function h = plotFullMesh(mesh_data, show_vertices, enable, usePastel)
            if nargin < 4, usePastel = false; end
            if nargin < 3, enable = true; end
            if nargin < 2, show_vertices = false; end
            if ~enable, h = []; return; end
            
            % --- Validaciones ---
            if ~isstruct(mesh_data) || ~isfield(mesh_data,'coords') || ~isfield(mesh_data,'connect')
                warning('mesh_data inválido. Se espera struct con coords y connect.');
                h = []; return;
            end
            coords = mesh_data.coords;
            if isempty(coords) || size(coords,2) ~= 2 || any(isnan(coords(:)) | isinf(coords(:)))
                warning('Coords inválidas'); h = []; return;
            end
            
            % Crear figura/axes (OpenGL recomendado)
            fig = figure('Renderer','opengl');
            ax = axes('Parent',fig);
            hold(ax,'on');
            
            % Para retorno
            h = struct();
            h.fig = fig;
            h.ax = ax;
            
            % --- Caso triángulos (mejor performance: un único patch) ---
            if ismatrix(mesh_data.connect) && size(mesh_data.connect,2) == 3
                faces = mesh_data.connect;
                verts = coords;
                
                % Colores pastel: generamos una paleta suave
                nFaces = size(faces,1);
                if usePastel
                    % Pastel: colores claros y desaturados
                    rng(0,'twister'); % reproducible; quita si quieres aleatorio
                    pastel = 0.75 + 0.25*rand(nFaces,3); % valores altos => claros
                    % opcionalmente reducir saturación aproximada
                else
                    pastel = repmat([1 1 1]*0.95, nFaces, 1); % casi blanco (solo borde visible)
                end
                
                % Crear patch con FaceColor por cara
                p = patch('Faces',faces,'Vertices',verts, ...
                    'FaceVertexCData',pastel, ...
                    'FaceColor','flat', ...    % flat permite un color por cara
                    'EdgeColor','k', ...
                    'LineWidth',0.4, ...
                    'Parent',ax, ...
                    'FaceLighting','none', ...
                    'EdgeAlpha',0.6);
                set(p,'HitTest','off','PickableParts','none');
                
                h.patch = p;
                
                % Mejora interacción: esconder edges mientras se hace zoom/pan
                storeEdge = get(p,'EdgeColor');
                zm = zoom(fig); pm = pan(fig);
                set(zm,'ActionPreCallback',@(obj,ev) set(p,'EdgeColor','none'));
                set(zm,'ActionPostCallback',@(obj,ev) set(p,'EdgeColor',storeEdge));
                set(pm,'ActionPreCallback',@(obj,ev) set(p,'EdgeColor','none'));
                set(pm,'ActionPostCallback',@(obj,ev) set(p,'EdgeColor',storeEdge));
                
                % --- Caso polígonos (cell array) ---
            elseif iscell(mesh_data.connect)
                % Montamos vectores con NaN separators y dibujamos en una sola llamada a line
                nPolys = numel(mesh_data.connect);
                x_all = [];
                y_all = [];
                faceColors = zeros(nPolys,3);
                validPolyIdx = false(1,nPolys);
                vc = 0;
                for i = 1:nPolys
                    idx = mesh_data.connect{i};
                    if numel(idx) < 3, continue; end
                    poly = coords(idx, :);
                    if ~isequal(poly(1,:), poly(end,:))
                        poly = [poly; poly(1,:)];
                    end
                    vc = vc + 1;
                    x_all = [x_all, poly(:,1)', NaN]; %#ok<AGROW>
                    y_all = [y_all, poly(:,2)', NaN]; %#ok<AGROW>
                    validPolyIdx(i) = true;
                    if usePastel
                        faceColors(i,:) = 0.80 + 0.2*rand(1,3);
                    else
                        faceColors(i,:) = [1 1 1]*0.95;
                    end
                end
                
                % Primero: crear patches (relleno pastel) si se pidió y hay polígonos
                if usePastel
                    % Si son muchos polígonos, advertir y crear igual (puede costar)
                    maxCreate = Inf; % pon Inf o un valor alto si quieres forzar
                    if vc > 1000
                        %warning('Muchas caras (%d). Relleno puede ser lento.', vc);
                    end
                    patches = gobjects(nPolys,1);
                    for i = 1:nPolys
                        if ~validPolyIdx(i), continue; end
                        idx = mesh_data.connect{i};
                        if numel(idx) < 3, continue; end
                        poly = coords(idx,:);
                        if ~isequal(poly(1,:), poly(end,:)), poly = [poly; poly(1,:)]; end
                        c = faceColors(i,:);
                        patches(i) = patch('XData',poly(:,1),'YData',poly(:,2),...
                            'FaceColor',c,'EdgeColor','none','FaceAlpha',0.9,'Parent',ax);
                        set(patches(i),'HitTest','off','PickableParts','none');
                    end
                    h.patches = patches;
                end
                
                % Luego: dibujar contornos en una sola llamada (siempre encima de fills)
                line_obj = line(x_all, y_all, 'Parent',ax, 'Color',[0.2 0.2 0.2], 'LineWidth',0.6);
                set(line_obj,'HitTest','off','PickableParts','none');
                h.line = line_obj;
                
                % Optimización interacción: ocultar patches/edges durante zoom/pan
                if usePastel
                    zm = zoom(fig); pm = pan(fig);
                    set(zm,'ActionPreCallback',@(obj,ev) set(line_obj,'Visible','off'));
                    set(zm,'ActionPostCallback',@(obj,ev) set(line_obj,'Visible','on'));
                    set(pm,'ActionPreCallback',@(obj,ev) set(line_obj,'Visible','off'));
                    set(pm,'ActionPostCallback',@(obj,ev) set(line_obj,'Visible','on'));
                end
            else
                warning('Formato de mesh_data.connect no soportado.');
                hold(ax,'off'); h = []; return;
            end
            
            % --- Vértices (opcional) ---
            if show_vertices
                s = scatter(ax, coords(:,1), coords(:,2), 8, 'MarkerFaceColor', 0.2*[1 1 1]+0.8, 'MarkerEdgeColor','none');
                set(s,'HitTest','off','PickableParts','none');
                h.scatter = s;
            end
            
            % Ajustes finales
            axis(ax,'equal');
            axis(ax,'tight');
            title(ax,'Malla Completa');
            hold(ax,'off');
            drawnow;
        end
    end
end