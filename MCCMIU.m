% Function: Obtaining the complete fuzzy sociometric and DMs'weights
clc,clear all;
format short g;
fprintf(2,'Calculation of the DMs weights£º\n')

%% 1.Input Data
data=[0    0.65 0    0.55 0.73 0    0.70 0    0.61;
      0.47 0    0.85 0    0    0.80 0    0.85    0;
      0.26 0.74 0    0.33 0    0.76 0.85 0.61 0.53;
      0.63 0    0.67 0    0.82 0    0    0    0.77;
      0    0    0    0.75 0    0    0    0.90    0;
      0.25 0.84 0.80 0    0    0    0.71 0       0;
      0.33 0    0.85 0.42 0.80 0.77 0    0.56    0;
      0    0.78 0    0    0.69 0    0.90 0    0.78;
      0.54 0.80 0.88 0.40 0    0.72 0    0.82    0];
[m,n]=size(data);
syms j1 j2 j3 j4 j5 j6 j7 j8 u1 u2 u3 u4 u5 u6 u7 u8
%% 2.Obtain complete fuzzy sociometric
Trust=[];  %Store the complete trust values of DMs
%2.1 Trust propagation in the single-path
for i=1:m 
    T=zeros(m,m-1);
    t=0;
    for j1=1:m
        t=1;
        x=[];
        if data(j1,i)>0
           x=[x,data(j1,i)];
           u1=x;
           T(j1,t)=EACTP(x); 
           for j2=1:m
               t=2;
               x=u1;
               if data(j2,j1)>0 && j2~=i && j2~=j1
                  x=[x,data(j2,j1)]; 
                  T(j2,t)=EACTP(x);
                  u2=x;
                  for j3=1:m
                      t=3;
                      x=u2;
                      if data(j3,j2)>0 && j3~=i && j3~=j1 && j3~=j2
                         x=[x,data(j3,j2)];
                         T(j3,t)=EACTP(x);
                         u3=x;
                         for j4=1:m
                             t=4;
                             x=u3;
                            if data(j4,j3)>0 && j4~=i && j4~=j1 && j4~=j2 && j4~=j3
                               x=[x,data(j4,j3)];
                               T(j4,t)=EACTP(x);
                               u4=x;
                               for j5=1:m
                                   t=5;
                                   x=u4;
                                   if data(j5,j4)>0 && j5~=i && j5~=j1 && j5~=j2 && j5~=j3 && j5~=j4
                                      x=[x,data(j5,j4)];
                                      T(j5,t)=EACTP(x);
                                      u5=x;
                                      for j6=1:m
                                          t=6;
                                          x=u5;
                                          if data(j6,j5)>0 && j6~=i && j6~=j1 && j6~=j2 && j6~=j3 && j6~=j4 && j6~=j5
                                             x=[x,data(j6,j5)];
                                             T(j6,t)=EACTP(x);
                                             u6=x;
                                             for j7=1:m
                                                 t=7;
                                                 x=u6;
                                                 if data(j7,j6)>0 && j7~=i && j7~=j1 && j7~=j2 && j7~=j3 && j7~=j4 && j7~=j5 && j7~=j6
                                                    x=[x,data(j7,j6)];
                                                    T(j7,t)=EACTP(x);
                                                    u7=x;
                                                    for j8=1:m
                                                        t=8;
                                                        x=u7;
                                                        if data(j8,j7)>0 && j8~=i && j8~=j1 && j8~=j2 && j8~=j3 && j8~=j4 && j8~=j5 && j8~=j6 && j8~=j7
                                                           x=[x,data(j8,j7)];
                                                           T(j8,t)=EACTP(x);
                                                           u8=x;
                                                        end
                                                    end
                                                 end
                                             end
                                          end
                                      end
                                   end
                                end
                            end
                         end
                      end
                  end
               end
           end
        end
    end
    
%2.2 Trust aggregation in multiple paths
    count=sum(T~=0,2);  
    Tr=zeros(m,1); 
    for i=1:m
        if T(i,1)>0 
           Tr(i,1)=T(i,1);
        else if count(i)==1  
                Tr(i,1)=sum(T(i,:),2);
             else if count(i)>1  
             %2.2.1 Determine the weight of each path
                  PL=[];
                  w=[];
                  for j=1:m-1
                      if T(i,j)>0
                         PL=[PL 1/j];
                      end   
                  end
                  for k=1:size(PL,2)
                      if k==1
                         v=(sum(PL(1:k))/sum(PL))^0.5; 
                         w=[w v];
                      else
                         v=(sum(PL(1:k))/sum(PL))^0.5-(sum(PL(1:k-1))/sum(PL))^0.5; 
                         w=[w v];
                      end
                  end
             %2.2.2 Determine the weight of each DM
                  P=T(i,:);
                  L=P(P~=0);
                  Tr(i,1)=w*L';
                 end
            end
        end
    end 
    Trust=[Trust Tr];
    
end
disp('The complete fuzzy sociometric:');
disp(Trust);

%% 3.Obtain the weights of DMs
W=[];
for i=1:m
    W(i)=sum(Trust(:,i),1)/sum(sum(Trust));
end
disp('Weights of DMs:');
disp(W);
