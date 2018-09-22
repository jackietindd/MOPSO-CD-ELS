% ----------------------------------------------------------------------- %
% Example of use of the funcion MOPSO.m, which performs a Multi-Objective %
% Particle Swarm Optimization (MOPSO), based on Coello2004.               %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez Cagigal                                      %
%   Date:    15/03/2017                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
%   References:                                                           %
%       Coello, C. A. C., Pulido, G. T., & Lechuga, M. S. (2004). Handling%
%       multiple objectives with particle swarm optimization. IEEE Tran-  %
%       sactions on evolutionary computation, 8(3), 256-279.              %
% ----------------------------------------------------------------------- %
%%==================================================
% example.m for MOPSOCD/CDELS
% Revised by shuxin ding
% version: 1.1
% last update: Sep 11th, 2018
% References:
% S. Ding, C. Chen, B. Xin, P M. Pardalos. A bi-objective load balancing
% model in a distributed simulation system using NSGA-II and MOPSO appro-
% aches. Applied Soft Computing, 2018, 63: 249-267.
%%==================================================


clear all; clc;

% Multi-objective function
%MultiObjFnc = 'Schaffer';
MultiObjFnc = 'Kursawe';
%MultiObjFnc = 'Poloni';
%MultiObjFnc = 'Viennet2';
%MultiObjFnc = 'Viennet3';
%MultiObjFnc = 'ZDT1';
%MultiObjFnc = 'ZDT2';
%MultiObjFnc = 'ZDT3';
%MultiObjFnc = 'ZDT6';

switch MultiObjFnc
    case 'Schaffer'         % Schaffer
        MultiObj.fun = @(x) [x(:).^2, (x(:)-2).^2];
        MultiObj.nVar = 1;
        MultiObj.var_min = -5;
        MultiObj.var_max = 5;
        load('ParetoFronts/Schaffer.mat');
        MultiObj.truePF = PF;
    case 'Kursawe'          % Kursawe 
        MultiObj.fun = @(x) [-10.*(exp(-0.2.*sqrt(x(:,1).^2+x(:,2).^2)) + exp(-0.2.*sqrt(x(:,2).^2+x(:,3).^2))), ...
                             sum(abs(x).^0.8 + 5.*sin(x.^3),2)];
        MultiObj.nVar = 3;
        MultiObj.var_min = -5.*ones(1,MultiObj.nVar);
        MultiObj.var_max = 5.*ones(1,MultiObj.nVar);
        load('ParetoFronts/Kursawe.mat');
        MultiObj.truePF = PF;
    case 'Poloni'           % Poloni's two-objective
        A1 = 0.5*sin(1)-2*cos(1)+sin(2)-1.5*cos(2);
        A2 = 1.5*sin(1)-cos(1)+2*sin(2)-0.5*cos(2);
        B1 = @(x,y) 0.5.*sin(x)-2.*cos(x)+sin(y)-1.5.*cos(y);
        B2 = @(x,y) 1.5.*sin(x)-cos(x)+2.*sin(y)-0.5.*cos(y);
        f1 = @(x,y) 1+(A1-B1(x,y)).^2+(A2-B2(x,y)).^2;
        f2 = @(x,y) (x+3).^2+(y+1).^2;
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)), f2(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min = -pi.*ones(1,MultiObj.nVar);
        MultiObj.var_max = pi.*ones(1,MultiObj.nVar);
    case 'Viennet2'         % Viennet2
        f1 = @(x,y) 0.5.*(x-2).^2+(1/13).*(y+1).^2+3;
        f2 = @(x,y) (1/36).*(x+y-3).^2+(1/8).*(-x+y+2).^2-17;
        f3 = @(x,y) (1/175).*(x+2.*y-1).^2+(1/17).*(2.*y-x).^2-13;
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)), f2(x(:,1),x(:,2)), f3(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min = [-4, -4];
        MultiObj.var_max = [4, 4];
        load('ParetoFronts/Viennet2.mat');
        MultiObj.truePF = PF;
    case 'Viennet3'         % Viennet3
        f1 = @(x,y) 0.5.*(x.^2+y.^2)+sin(x.^2+y.^2);
        f2 = @(x,y) (1/8).*(3.*x-2.*y+4).^2 + (1/27).*(x-y+1).^2 +15;
        f3 = @(x,y) (1./(x.^2+y.^2+1))-1.1.*exp(-(x.^2+y.^2));
        MultiObj.fun = @(x) [f1(x(:,1),x(:,2)), f2(x(:,1),x(:,2)), f3(x(:,1),x(:,2))];
        MultiObj.nVar = 2;
        MultiObj.var_min = [-3, -10];
        MultiObj.var_max = [10, 3];
        load('ParetoFronts/Viennet3.mat');
        MultiObj.truePF = PF;
    case 'ZDT1'             % ZDT1 (convex)
        g = @(x) 1+9.*sum(x(:,2:end),2)./(size(x,2)-1);
        MultiObj.fun = @(x) [x(:,1), g(x).*(1-sqrt(x(:,1)./g(x)))];
        MultiObj.nVar = 30; 
        MultiObj.var_min = zeros(1,MultiObj.nVar);
        MultiObj.var_max = ones(1,MultiObj.nVar);
        load('ParetoFronts/ZDT1.mat');
        MultiObj.truePF = PF;
    case 'ZDT2'             % ZDT2 (non-convex)
        f = @(x) x(:,1);
        g = @(x) 1+9.*sum(x(:,2:end),2)./(size(x,2)-1);
        h = @(x) 1-(f(x)./g(x)).^2;
        MultiObj.fun = @(x) [f(x), g(x).*h(x)];
        MultiObj.nVar = 30; 
        MultiObj.var_min = zeros(1,MultiObj.nVar);
        MultiObj.var_max = ones(1,MultiObj.nVar);
        load('ParetoFronts/ZDT2.mat');
        MultiObj.truePF = PF;
    case 'ZDT3'             % ZDT3 (discrete)
        f = @(x) x(:,1);
        g  = @(x) 1+(9/size(x,2)-1).*sum(x(:,2:end),2);
        h  = @(x) 1 - sqrt(f(x)./g(x)) - (f(x)./g(x)).*sin(10.*pi.*f(x));
        MultiObj.fun = @(x) [f(x), g(x).*h(x)];
        MultiObj.nVar = 30;
        MultiObj.var_min = 0.*ones(1,MultiObj.nVar);
        MultiObj.var_max = 1.*ones(1,MultiObj.nVar);
        load('ParetoFronts/ZDT3.mat');
        MultiObj.truePF = PF;
    case 'ZDT6'             % ZDT6 (non-uniform)
        f = @(x) 1 - exp(-4.*x(:,1)).*sin(6.*pi.*x(:,1));
        g = @(x) 1 + 9.*(sum(x(:,2:end),2)./(size(x,2)-1)).^0.25;
        h = @(x) 1 - (f(x)./g(x)).^2;
        MultiObj.fun = @(x) [f(x), g(x).*h(x)];
        MultiObj.nVar = 10;
        MultiObj.var_min = 0.*ones(1,MultiObj.nVar);
        MultiObj.var_max = 1.*ones(1,MultiObj.nVar);
        load('ParetoFronts/ZDT6.mat');
        MultiObj.truePF = PF;
end

% Parameters
params.Np = 200;        % Population size
params.Nr = 200;        % Repository size
params.maxgen = 500;    % Maximum number of generations
params.W = 0.4;         % Inertia weight
params.C1 = 2;          % Individual confidence factor
params.C2 = 2;          % Swarm confidence factor
params.maxvel = 5;      % Maxmium vel in percentage

% MOPSO
REP = MOPSOCDELS(params,MultiObj);

% Display info
display('Repository fitness values are stored in REP.pos_fit');
display('Repository particles positions are store in REP.pos');


    