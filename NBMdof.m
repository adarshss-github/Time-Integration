function [Xd,Xv,Xa,t] = NBMdof(alpha,beta,M,C,K,F,fs,Xd0,Xv0,tf)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[Xd,Xv,Xa,t] = NBMdof(alpha,beta,M,C,K,F,fs,Xd0,Xv0,tf)
%
%Adarsh S, Ph.D. Candidate, IIT Kanpur
%
%Description:
%-----------
%This function calculates the displacement, velocity and acceleration
%vectors of an LTI mdof system, given the excitation, using the generalized
%Newmark Beta time stepping algorithm
%
%Input arguments:
%---------------
%1)alpha: Newmark beta constant
%2)beta: Newmark beta constant
%3)M:Mass matrix of the mdof system
%4)C: Damping matrix of the mdof system
%5)K:Stiffness matrix of the mdof system
%6)F: Matrix containing the forces at each dof, arranged in row wise manner
%7)fs: Sampling frequency
%8)Xd0: Vector of initial displacements at each dof
%9)Xv: Vector of initial velocities at each dof
%10)tf: Final time upto which you want to find the response in sec.
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
%
%Ex:  [Xd,Xv,Xa,t] = NBMdof(1/2,1/6,M,C,K,F,200,Xd0,Xv0,40) ;
%
%Note:
%---- 
%1) If alpha = 1/2 and beta = 0, then this will be the Central difference algorithm
%2) If alpha = 1/2 and beta = 1/4, then this will be the Constant acceleration algorithm
%3) If alpha = 1/2 and beta = 1/6, then this will be the Linear acceleration algorithm
%4) This code can be used for sdof systems also, where M,C,K,Xd0,Xv0 reduce
%to scalars and F will become a row vector
%5) You may also use it for structures subjected to ground motion,
%   in this case F = -M*{influence vector}*uga
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; November, 2017 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================


[~,D] = eig(K,M) ;
f = eigtoHz(D) ;
T = 1/max(f) ;
dt = 1/fs ;
clear D ;

%Stability Check
if alpha>= 1/2 && beta< alpha/2 && dt/T > ( 1/( 2*pi*sqrt( (alpha/2) - beta ) ) ) 
    
    error('Newmark beta algorithm wont be stable, consider increasing fs ') ;
    
end


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

D = M + ((alpha)*dt)*C + ((beta)*dt^2)*K ;
iD = inv(D) ;
%L = chol(D,'lower') ; If you want to use Cholesky factorization
clear D ;
P = zeros(mF,1) ;

for i = 1:1:N-1
   
    P = F(:,i+1) - C*( Xv(:,i) + ( 1-(alpha) )*dt*Xa(:,i) ) - K*( Xd(:,i) + dt*Xv(:,i) + ( (1/2) - beta )*dt^2*Xa(:,i) ) ;
    %Xa(:,i+1) = L'\(L\P) ; If you want to use Cholesky factorization
    Xa(:,i+1) = iD*P ;
    Xv(:,i+1) = Xv(:,i) + (1-alpha)*dt*Xa(:,i) + alpha*dt*Xa(:,i+1) ;
    Xd(:,i+1) = Xd(:,i) + dt*Xv(:,i) + ((1/2) - beta )*dt^2*Xa(:,i) + beta*dt^2*Xa(:,i+1) ;
    
end    
 
clear alpha beta M C K F fs Xd0 Xv0 tf dt f T P iD
    
end

function [f] = eigtoHz(D)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
% [f] = eigtoHz(D)
%
%
%Adarsh S, Ph.D. Candidate, IIT Kanpur
%
%Description:
%-----------
%This function takes the digonal matrix of eigen values and returns the natural frequencies in Hz.
%
%Input arguments:
%---------------
%1)D: Diagonal matrix of eigen values
%
%Output arguments:
%----------------
%1)f: Natural frequencies in Hz.
%
%Ex:  
%[~,D] = eig(K,M) ;
%f = eigtoHz(D) ;
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; November, 2017 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

w2 = diag(D) ;
 w = sqrt(w2) ;
 f = w/(2*pi) ; % Natural frequencies in Hz
 clear w2 w ;

end
    