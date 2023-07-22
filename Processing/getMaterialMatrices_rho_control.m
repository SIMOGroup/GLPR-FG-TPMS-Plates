function [Ab,Bb,Db,Eb,Fb,Hb,Ds,I]=getMaterialMatrices_rho_control(material,h,imix)
%%
global rho_distribution M_mean RD_max RD_0
global GPL_distribution l_gpl w_gpl t_gpl E_gpl nu_gpl rho_gpl lambda_gpl
global porous_type

% ===== Material matrix =====
Ab=zeros(3,3); Bb=Ab; Db=Ab; Eb=Ab; Fb=Ab; Hb=Ab;
Ds=zeros(2,2);I=zeros(1,6);

% % ============ material property============================
% get integration points
%[Wg,Qg] = quadrature(8,'GAUSS',1)

x1 = 1;
x2 = -1 ;
[Qg,Wg]=gauleg(x1,x2,20) ;
zu = h/2;              % upper
zl = -h/2;             % lower

% Young's modulus, poisson's ratio
E_m = material(1);
nu_m = material(2);
G_m = E_m/(2*(1+nu_m));
rho_m = material(3);

%% Define parameter of TPMS
switch porous_type
    case {1, 2, 3}
        if porous_type == 1 % Primitive
            k_e = 0.25; C_1e = 0.317; n_1e = 1.264; n_2e = 2.006;
            k_g = 0.25; C_1g = 0.705; n_1g = 1.189; n_2g = 1.715;
            k_nu = 0.55; a_1 = 0.314; b_1 = -1.004; a_2 = 0.152;
        elseif porous_type == 2 % Gyroid
            k_e = 0.45; C_1e = 0.596; n_1e = 1.467; n_2e = 2.351;
            k_g = 0.45; C_1g = 0.777; n_1g = 1.544; n_2g = 1.982;
            k_nu = 0.50; a_1 = 0.192; b_1 = -1.349; a_2 = 0.402;
        elseif porous_type == 3 % IWP
            k_e = 0.35; C_1e = 0.597; n_1e = 1.225; n_2e = 1.782;
            k_g = 0.35; C_1g = 0.529; n_1g = 1.287; n_2g = 2.188;
            k_nu = 0.13; a_1 = 2.597; b_1 = -0.157; a_2 = 0.201;
        end
        C_2e = (C_1e*k_e^(n_1e) - 1)/(k_e^(n_2e) - 1); C_3e = 1 - C_2e;
        C_2g = (C_1g*k_g^(n_1g) - 1)/(k_g^(n_2g) - 1); C_3g = 1 - C_2g;
        d_1 = 0.3 - a_1*exp(b_1*k_nu); b_2 = - a_2*(k_nu + 1); d_2 = 0.3 - a_2*(1)^2 - b_2(1);
end
            
%% Calculate n_{A}, n_{B}, psi_{rho} of porous ditribution 1, 2, 3
switch rho_distribution
    case 1
        n_A = (RD_max - M_mean) / (M_mean - RD_max*(1-RD_0));
    case 2
        n_B = 1;
        RD_min = RD_max*(1- RD_0);
        try
            csv_data = csvread(['rho_framework/0_Max_', char(string(M_mean*100)), '.csv']);
            index_0 = find(abs(csv_data(:,1) - RD_0) <= 1e-3, 1, 'last');
            index_max = find(abs(csv_data(1,:) - RD_max) <= 1e-3, 1, 'last');
            n_B = csv_data(index_0, index_max);
            if isnan(n_B)
                disp("Error in finding n_B")
                pause
            end
        catch
            disp('Error in finding .csv');
            pause
        end
    case 3
        psi_rho = (M_mean - RD_max*(1-RD_0)) / (RD_max*RD_0);
end

%% Calculate Si of GPL distributions
Vtotal_GPL = lambda_gpl*rho_m / (lambda_gpl*rho_m + rho_gpl - lambda_gpl*rho_gpl);
A1 = 0; A2 = 0;
for igp = 1:size(Wg,1)
    pt = Qg(igp,:) ;
    
    % map the point to global
    gpt = (zu+zl)/2 + (zu-zl)/2*pt ;                 % value of z-direction
    wt = (zu-zl)/2*Wg(igp) ;
    
    switch rho_distribution
        case 1
            psi_z = (gpt/h + 1/2)^n_A;
        case 2
            psi_z = (1 - cos(pi*gpt/h))^n_B;
        case 3
            psi_z = psi_rho;
    end
    
    RD_z = RD_max * (1-RD_0 + RD_0*psi_z);

    switch GPL_distribution
        case 1
            Vgpl_z = (1 - cos(pi*gpt/h));
        case 2
            Vgpl_z = (1 - cos(pi*gpt/2/h + pi/4));
        case 3
            Vgpl_z = 1;
    end
    
    A1 = A1 + RD_z*wt;
    A2 = A2 + (Vgpl_z)*(RD_z)*wt;
end
S = Vtotal_GPL*A1/A2;

%% Calculate material matrices
xi_l = 2*l_gpl/t_gpl; xi_w = 2*w_gpl/t_gpl;
eta_l = ((E_gpl/E_m) - 1)/((E_gpl/E_m) + xi_l); eta_w = ((E_gpl/E_m) - 1)/((E_gpl/E_m) + xi_w);
for igp = 1:size(Wg,1)
    pt = Qg(igp,:) ;
    
    % map the point to global
    gpt = (zu+zl)/2 + (zu-zl)/2*pt ;                 % value of z-direction
    wt = (zu-zl)/2*Wg(igp) ;
    
    % GPLR base material
    switch GPL_distribution
        case 1
            Vgpl_z = S * (1 - cos(pi*gpt/h));
        case 2
            Vgpl_z = S * (1 - cos(pi*gpt/2/h + pi/4));
        case 3
            Vgpl_z = S;
    end
    
    V_m = 1 - Vgpl_z;
    E_1 = 3/8*(1 + xi_l*eta_l*Vgpl_z)/(1-eta_l*Vgpl_z)*E_m + 5/8*(1 + xi_w*eta_w*Vgpl_z)/(1-eta_w*Vgpl_z)*E_m;
    nu_1 = nu_gpl*Vgpl_z + nu_m*V_m;
    G_1 = E_1 / (2*(1+nu_1));
    rho_1 = rho_gpl*Vgpl_z + rho_m*V_m;
    
    % Including the porosity
    switch rho_distribution
        case 1
            psi_z = (gpt/h + 1/2)^n_A;
        case 2
            psi_z = (1 - cos(pi*gpt/h))^n_B;
        case 3
            psi_z = psi_rho;
    end
    
    % 1: Primitive, 2: Gyroid, 3: IWP, 4: Closed-cell, 5: Open-cell
    RD_z = RD_max * (1-RD_0 + RD_0*psi_z);
    switch porous_type
        case {1, 2, 3}
            e_z = (RD_z <= k_e) * (C_1e*RD_z^(n_1e)) + ...
                  (RD_z >  k_e) * (C_2e*RD_z.^(n_2e) + C_3e);
            g_z = (RD_z <= k_g) * (C_1g*RD_z^(n_1g)) + ...
                  (RD_z >  k_g) * (C_2g*RD_z.^(n_2g) + C_3g);
            nu_z =(RD_z <= k_nu) * (a_1.*exp(b_1*RD_z) + d_1) + ...
                  (RD_z >  k_nu) * (a_2*RD_z^(2) + b_2*RD_z + d_2);
            
            E_z = E_1*e_z;
            G_z = G_1*g_z;
%             rho_z = rho*RD;
        case 4
            E_z = E_1*((RD_z + 0.121)/1.121)^2.3;
            nu_z = 0.221*(1-RD_z) + nu_1*(0.342*(1-RD_z)^(2) - 1.21*(1-RD_z) + 1);
            G_z = E_z / (2*(1+nu_z));
        case 5
            E_z = E_1*(RD_z)^2;
            nu_z = 1/3;
            G_z = E_z /(2*(1+nu_z));
        case 6
            E_z = E_1*(RD_z)^2;
            nu_z = 0.3;
            G_z = E_z /(2*(1+nu_z));
    end
    rho_z = rho_1 * RD_z;
    
    % % ============= Bending-extension stiffness matrix ===================
    Qb = [E_z/(1-nu_z^2) E_z*nu_z/(1-nu_z^2) 0;
          E_z*nu_z/(1-nu_z^2) E_z/(1-nu_z^2) 0;
          0 0 G_z];
    % shear stiffness
    Qs =G_z*[1 0;
             0 1];
    % choos model for displacement field
    switch imix
        case 1
            %E. Reissner, The effect of transverse shear deformation on the bending of elastic plates
            gg = 0;% 3DOF with no higher-order term %
            dgg= 0;
        case 2
            %J. N. Reddy, A Simple Higher-Order Theory for Laminated Composite Plates,
            gg = gpt-4*gpt^3/(3*h^2);% Reddy
            dgg= 1-4*gpt^2/h^2;
        case 3
            %R. P. Shimpi, Refined Plate Theory and Its Variants
            gg = gpt-5*gpt^3/(3*h^2)+gpt/4;%Shimpi
            dgg= 1-5*gpt^2/h^2+1/4;
        case 4
            %H. Nguyen-Xuan, Isogeometric finite element analysis of composite sandwich plates using a higher order shear deformation theory
            gg = gpt-1/8*gpt-2/h^2*gpt^3+2/h^4*gpt^5;% H. Nguyen-Xuan
            dgg= 1-1/8-6/h^2*gpt^2+10/h^4*gpt^4;
        case 5
            % H. X. Nguyen, T. N. Nguyen, A refined quasi-3D isogeometric analysis for functionally graded microplates based on the modified couple stress theory
            gg = gpt-9*gpt+(10/h^2)*gpt^3+6/(5*h^4)*gpt^5+8/(7*h^6)*gpt^7;%Hoang X. Nguyen
            dgg= 1-9+(10/h^2)*3*gpt^2+6/(5*h^4)*5*gpt^4+8/(7*h^6)*7*gpt^6;
        case 6
            %T. N. Nguyen, On the general framework of high order shear deformation theories for laminated composite plate structures: A novel unified approach
            gg = gpt - (17*gpt^3)/(10*h^2) + (22*gpt^5)/(25*h^4);%Tuan N. Nguyen
            dgg=  1- (17*3*gpt^2)/(10*h^2) + (22*5*gpt^4)/(25*h^4);
        case 7
            %C. H. Thai, Generalized shear deformation theory for functionally graded isotropic and sandwich plates based on isogeometric approach
            gg = atan(sin(pi/h*gpt));%Chien H. Thai (Computers & Structures)
            dgg= pi/h*cos(pi/h*gpt)/(1+(sin(pi/h*gpt))^2);
        otherwise
            disp('**************************************************');
            disp('Do not appropriate models ');
            disp('**************************************************');
            pause
    end
    
    Ab = Ab + Qb*wt ;        Bb = Bb + Qb*gpt*wt ;     Db = Db + Qb*gpt^2*wt ;
    Eb = Eb + Qb*gg*wt ;     Fb = Fb + Qb*gg*gpt*wt ;  Hb = Hb + Qb*gg^2*wt ;
    Ds = Ds + Qs*(dgg)^2*wt;         
    % ================= Inertia terms matrix =======================
    I = I + (rho_z*wt).* [1 gpt gpt^2 gg gg*gpt gg^2];
end