% **********************************************************************
% *                                                                    *
% *     This program determines the response of a multistory single    *
% *     bay moment frame building (w/ infinitely rigid beams) under    *
% *     a given ground motion.                                         *
% *                     |                 |                            *
% *                     O======( mi )=====O                            *
% *                     |                 |                            *
% *                k_i  |                 |  k_i                       *
% *                     |                 |                            *
% *                     O=====( m_i-1)====O                            *
% *                     |                 |                            *
% *                                                                    *
% *                                                                    *
% *                     |                 |                            *
% *             node 3  O======( m1 )=====O  node 4                    *
% *                     |       ele3      |                            *
% *                k_1  | ele1        le2 |  k_1                       *
% *                     |                 |                            *
% *            node1  ~~~~~             ~~~~~  node 2                  *
% *               <---------->      <---------->                       *
% *                   ag(t)             ag(t)                          *
% *                                                                    *
% *-----INPUT:                                                         *
% *                                                                    *
% *      <dt>    : Time increment                                      *
% *      <ag>    : nsteps x 1 vector of ground motion accelerations    *
% *      <kfl>   : ndof x 1 vector of floor (bending) stiffnesses      *
% *      <mfl>   : ndof x 1 vector of floor masses                     *
% *      <kdamp> : Stiffness proportional damping coefficient          *
% *      <mdamp> : Mass proportional damping coefficient               *
% *                                                                    *
% *-----OUTPUT:                                                        *
% *                                                                    *
% *      <u,v,a> : ndof x nsteps matrices of floor displacements,      *
% *                velocities and accelerations                        *
% *      <M,C,K> : ndof x ndof matrices of mass, damping and stiffness * 
% *      <w>     : ndof x 1 vector of modal frequencies                * 
% *      <xi>    : ndof x 1 vector of modal damping ratios             * 
% *      <Phi>   : ndof x ndof matrix of mode shapes                   * 
% *                                                                    *
% *-----GLOSSARY:                                                      *
% *                                                                    *
% *      <ndof>  : Number of degrees-of-freedom (floors)               *
% *                                                                    *
% **********************************************************************
% *                                                                    *
% *   MomentFrameBuilding                                              *
% *                                                                    *
% *                                                                    *
% *                                                                    *
% **********************************************************************

close all
clear all
clc
  
% Give the project a title
  ptitle =  'Six-Story Moment Frame Building'

% El Centro NS ground motion
  
  load ElCentroNS.txt 
  ag90 = column (ElCentroNS);
     
% Set units for ground motion
  ag = 386.4 * ag90;  % north-south acceleration                   
  clear ag90 ; clear El_Centro_Chopra_194x8.txt; % Remove raw data from memory
  dt = 1/50;   
  time = dt:dt:dt*length(ag);  % Time Array

  % Number of stories
  nStories = 6;
  
  % Floor heights
  h_typ  = 10.5*12; 
  hfl = h_typ*ones(1, nStories);
  hfl(1) = 15*12; % Update first story to be 15 ft
  Htot = sum(hfl);

  % Create vector for story heights for y axis values
  modeheights = zeros(numel(hfl)+1,1);
  for i = 1:numel(hfl)+1
      if i == 1
          modeheights(i) = 0;
      elseif i == 2
          modeheights(i) = hfl(1);
      else
          modeheights(i) = hfl(i-1) + modeheights(i-1);
      end
  end
  modeheights = modeheights/12; % convert from in to feet
  
  nmodes = length(hfl);
  M    = zeros(nmodes);   
  C    = zeros(nmodes);  
  K    = zeros(nmodes);   
  w    = zeros(nmodes,1); 
  xi   = zeros(nmodes,1);

  % Generate Mass and stiffness matricies with functions
  M = Mc1NS();
  K = Kc1NS();

  disp('The mass matrix is:');
  disp(M);
  disp('The stiffness matrix is:');
  disp(K);

  % Modal analysis, obtain mode shapes and modal periods
  [eigVectors,eigValues]=eig(K,M);
  % eigValues is a diagonal matrix of eigenvalues
  [numOfModes,~]=size(K);
    
  [eigValues10x1,order]=sort(diag(eigValues),'ascend');
  Phin=eigVectors(:,order);

  fprintf('The Eigenvalues (sorted smallest to largest):\r');
  eigValues10x1
  fprintf('The Eigenvectors (corresponding w/ sorted E Values):\r');
  Phin
    
  plotModeShapes2(Phin,(modeheights),1,0);

  % Determine modal frequencies
  w=zeros(1,numOfModes);
  for c=1:numOfModes
      w(1,c)=sqrt(eigValues10x1(c,1));
  end
  T=(2*pi)./w;

  % Plot of the mode shapes
  color = hsv(nmodes);  level = [0 1:nmodes];  figure(2);
  for i = 1:nmodes;
      plot([0; Phin(:,i)],level,'color',color(i,:));
      hold on;
  end
  grid on
  clear leg
  for i=1:nmodes
      T_local = T(i);
      T_local = num2str(T_local);
      str = "Mode " + i + " (T = " + T_local + " s)";
      leg(i,:) = str;
  end  
  legend(leg)
  title('Mass Normalized Mode Shapes')
  xlabel('Displacement');  ylabel('Level');  hold off

  % Rayleigh Damping coefficients and matrix
  dratio = 0.05; % Set damping to 5% to allow comparison with Chopra
  a1 = 2*dratio/(w(1)+w(4)); % 5% @ mode 1&4
  a0 = a1*w(1)*w(4);
  xi = zeros(nmodes, 1);

  for i=1:nmodes
      xi(i)=a0/(2*w(i))+a1*w(i)/2;
  end
 
  C = a0*M+a1*K;
  CN = Phin'*C*Phin; % C normalized by phi
  cn = diag(CN); % Diagonal values of CN
  
  fprintf('The Rayleigh Damping coefficients (Damping Ratios) for all modes:\r'); xi
  fprintf('The Rayleigh Damping coefficients a0 and a1:\r'); a0, a1
  fprintf('The Rayleigh Damping Matrix:\r'); C
  fprintf('The modal damping matrix PhiT*C*Phi:\r'); CN
  fprintf('The diagonal modal damping values from PhiT*C*Phi:\r'); cn