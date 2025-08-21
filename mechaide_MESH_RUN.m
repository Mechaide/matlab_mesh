% =============================================================================
% MECHAIDE Matlab
% Version 1.0 [21-08-2025].
% -----------------------------------------------------------------------------
% File: mechaide_MESH_RUN.m v1.0 [21-08-2025].
% Description: Programa para cargar
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

clear;
clc;

% Parámetros
dxf_file = 'ejemplos/test1.dxf'; % Ruta al archivo DXF
d = 0.5; % Separación entre vértices
mesh_type = 'dpolylla'; % Tipo de malla: 'triangle' o 'dpolylla'
decimales = 2; % Decimales para redondeo de coordenadas
iterations = 3; % Iteraciones de relajación de Lloyd
plot_contour = true; % Graficar contorno
plot_base_mesh = true; % Graficar malla base
plot_frontier = true; % Graficar fronteras
plot_contour_vertices = true; % Graficar vértices del contorno
plot_mesh = true; % Graficar malla final

% Generar malla
mesh_data = mechaide_MESH.generateMeshFromDXF(dxf_file, d, mesh_type, ...
    decimales, iterations, plot_contour, plot_base_mesh, plot_frontier, plot_mesh);

