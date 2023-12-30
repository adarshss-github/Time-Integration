function [Xd,Xv,Xa,t] = HHTalpha(alpha,M,C,K,F,fs,Xd0,Xv0,tf)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[Xd,Xv,Xa,t] = HHTalpha(alpha,M,C,K,F,fs,Xd0,Xv0,tf) ;
%
%Adarsh S, Ph.D. Candidate, IIT Kanpur
%
%Function description:
%---------------------
%This function calculates the displacement, velocity and acceleration
%vectors of an LTI mdof system, given the excitation, using the
%Hilber-Hughes-Taylor alpha method a.k.a HHT-alpha method. This is useful
%when we want to damp out contributions of high frequency modes
%
%Input arguments:
%---------------
%1)alpha:Decides how much of the higher modes need to be damped out (Should be between 0 and 1/3, higher values imply more numerical damping) 
%2)M:Mass matrix of the mdof system
%3)C: Damping matrix of the mdof system
%4)K:Stiffness matrix of the mdof system
%5)F: Matrix containing the forces at each dof, arranged in row wise manner
%5)fs: Sampling frequency
%7)Xd0: Vector of initial displacements at each dof
%8)Xv: Vector of initial velocities at each dof
%9)tf: Final time upto which you want to find the response in sec.
%
%Output arguments:
%----------------
%1)Xd: Matrix containing the relative displacements at each dof, arranged in row
%wise manner
%2)Xv: Matrix containing the relative velocities at each dof, arranged in row
%wise manner
%3)Xa: Matrix containing the relative accelerations at each dof, arranged in row
%wise manner
%4)t: Time axis in sec.
%
%Ex:  [Xd,Xv,Xa,t] = HHTalpha(0.10,M,zeros(2,2),K,F,100,zeros(2,1),zeros(2,1),200) ;
%
%Note:
%---- 
%1) This algorithm is designed to be unconditionally stable with atleast
%second order accuracy; the Newmark beta linear acc. had second order
%accuracy but was conditionally stable and did not have numerical damping
%2) If aplha = 0, it reduces to Newmark beta algorithm with constant
%acceleration; increase the value of alpha for more numerical damping
%3) This code can be used for sdof systems also, where M,C,K,Xd0,Xv0 reduce
%to scalars and F will become a row vector
%4) You may also use it for structures subjected to ground motion,
%   in this case F = -M*{influence vector}*uga
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; November, 2017 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

if alpha < 0 || alpha > (1/3) 
    
    error('The value of alpha should lie between 0 and 1/3 ') ;
    
end

beta = ((1+ alpha)^2)/4 ;
gamma = (1/2) + alpha ;

dt = 1/fs ;
t = 0:dt:tf ;
N = length(t) ;
[mF,nF] = size(F) ;

if N>nF
    
   F = [F zeros(mF,N-nF)] ;
  [mF,nF] = size(F) ;
       
end

Xd = zeros(mF,N) ;
Xv = zeros(mF,N) ;
Xa = zeros(mF,N) ;

Xd(:,1) = Xd0 ;
Xv(:,1 ) = Xv0 ;
Xa(:,1) = inv(M)*( F(:,1) - C*Xv(:,1) - K*Xd(:,1) ) ;

dA = (M +dt*(1-alpha)*gamma*C + (dt^2)*(1-alpha)*beta*K) ;
dB = (dt*(1-alpha)*(1-gamma)*C + (dt^2)*(1-alpha)*(0.5-beta)*K) ;
dC = (C + dt*(1-alpha)*K) ;
invdA = inv(dA) ;
clear A ;

for i = 1:1:N-1
   
    Xa(:,i+1) = invdA*( -dB*Xa(:,i) -dC*Xv(:,i) -K*Xd(:,i) + (1-alpha)*F(:,i+1) + alpha*F(:,i) ) ;
    Xv(:,i+1) = Xv(:,i) + (1-gamma)*dt*Xa(:,i) + gamma*dt*Xa(:,i+1) ;
    Xd(:,i+1) = Xd(:,i) + dt*Xv(:,i) + ((1/2) - beta )*dt^2*Xa(:,i) + beta*dt^2*Xa(:,i+1) ;
    
end    

end