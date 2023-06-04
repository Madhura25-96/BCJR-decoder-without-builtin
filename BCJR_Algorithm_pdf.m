%%BCJR Algorithm 
clear;
trellis = poly2trellis(3,[6 7]);
output_num = log2(trellis.numOutputSymbols);
input_num = log2(trellis.numInputSymbols);
%%Coding Rate, R = input_num / output_num;
Rate = input_num/output_num;
mLen = log2(trellis.numStates);
k = mLen*output_num;
snr = 3;
%Output of Sample Data 101; with given trellis (3,[6 7]);
y = [1 1 1 1 1 0];
length_encMsg = length(y);
expMsg = length_encMsg*Rate;
decoded_msg = zeros(1,expMsg);
States = trellis.numStates;
%%Step1 compute gamma:- Branch Metric, sigma is estimated
sigma = 10^(-snr/10);
awgn_output = awgn(y,snr);
%%disp('AWGN Output');
%%disp(awgn_output);
Gamma = zeros(States,States,expMsg);
disp('TRELLLIS Outputs');
disp(trellis.outputs);
disp('Trellis Next States');
disp(trellis.nextStates);
%%Z = No of Stages.
z = 1;
prev_st = 1;
for present_st = 1:States
    %%If there is any branch
    [data2,index] =  ismember(present_st-1, trellis.nextStates(prev_st,:));
    if data2 == 1
       e = trellis.outputs(prev_st,index);
       %%this is done because we need 00 and 01 instead of and 1
        if (e == 0)
          e = [0 0];
        elseif (e == 1)
          e = [0 1];
        else
          e = oct2poly(e);
        end
          disp('z');
          disp(z);
          disp('prev_st');
          disp(prev_st);
          disp('present_st');
          disp(present_st);
          Gamma(prev_st,present_st,1)= 1*exp(2.5*(e(1)*awgn_output(2*z-1)+e(2)*awgn_output(2*z)));
          disp('Gamma');
          disp(Gamma);
    end
end
%%For prev_st = 2 and 3 we dont want to calculate
for z = 2:expMsg
            for prev_st = 1:States
                %%If for a node there is no branch, it means its sum from previous stage is 0, this way its gamma is zero
                %%gamma is zero even for present stage
                %%Example in stage 1 node 2 has no branch, so sum of all previous states[(1,2,1),(2,2,1),(3,2,1),(4,2,1)] = 0
                %%Gamma [(2,1,2),(2,2,2),(2,3,2),(2,4,2)] == 0
                if(sum(Gamma(:,prev_st,z-1)) == 0 )
                    Gamma(prev_st,:,z) == 0;
                else
                  for present_st = 1:States
                    [data2,index] =  ismember(present_st-1, trellis.nextStates(prev_st,:));
                    if data2 == 1
                        e = trellis.outputs(prev_st,index);
                 %%this is done because we need 00 and 01 instead of and 1
                        if (e == 0)
                            e = [0 0];
                        elseif (e == 1)
                            e = [0 1];
                        else
                            e = oct2poly(e);
                        end
                        disp('z');
                        disp(z);
                        disp('prev_st');
                        disp(prev_st);
                        disp('present_st');
                        disp(present_st);
                        Gamma(prev_st,present_st,z) = 1*exp(2.5*(e(1)*awgn_output(2*z-1)+e(2)*awgn_output(2*z)));
                        disp('Gamma');
                        disp(Gamma);
                  end
                end
            end 
end
end
%%Step 2:-  Alpha recursions
alpha = zeros(States,expMsg+1);
disp('Estimated Msg');
disp(expMsg);
disp('States');
disp(States);
alpha(:,1) = 0;
alpha(1,1) = 1;
for index = 2:expMsg+1
     %%alpha what matters is the prev state
     for present_st = 1:States
         for prev_st = 1:States
             %%alpha recursion
             alpha(present_st,index) = alpha(present_st,index) + Gamma(prev_st,present_st,index-1)*alpha(prev_st,index-1);
          end
     end
     alpha(:,index) = alpha(:,index)/sum(alpha(:,index)); % Normalization
 end
 %Step 3 :- ******beta recursions******
 beta = zeros(States,expMsg+1);
 beta(:,expMsg+1) = alpha(:,1);
 %%since we want backward for beta
 for inde=expMsg+1:-1:2
       for present_st = 1:States
           for prev_st = 1:States
                %Beta recursion
               beta(present_st,inde-1) = beta(present_st,inde-1) + Gamma(prev_st,present_st,inde-1)*beta(prev_st,inde);
           end
       end
      beta(:,inde-1) = beta(:,inde-1)/sum((beta(:,inde-1))) ; % Normalization
 end
 disp('Alphaa');
 disp(alpha);
 disp('Beta');
 disp(beta);
 %Step 4 :- Log likelihood is based on log( input 0/ input = 1)
 for q = 1:expMsg
       cnt1 = 0;
       cnt2 =0;
       JointProb_input1 = 0;
       JointProb_input2 = 0;
       for present_st = 1:States
              for prev_st = 1:States
                  [data2,index2] = ismember(prev_st-1,trellis.nextStates(present_st,:));
                 %%disp('data2');
                 %%disp(data2);
                 %%disp('index2');
                 %%disp(index2);
                  if (data2==1 && index2==1) %input = 0
                       JointProb_input1= JointProb_input1 + alpha(present_st,q)*Gamma(present_st,prev_st,q)*beta(prev_st,q+1);
                       disp('Joint Probability input1');
                       disp(JointProb_input1);
                       cnt1 = cnt1+1;
                  elseif (data2==1 && index2==2) %input = 1
                       JointProb_input2= JointProb_input2 + alpha(present_st,q)*Gamma(present_st,prev_st,q)*beta(prev_st,q+1);
                       disp('Joint Probability input2');
                       disp(JointProb_input2);
                       cnt2 = cnt2+1;
                  end
          end
          JointProb_input1 = JointProb_input1/(JointProb_input1 + JointProb_input2) ; % Normalization       
          JointProb_input2 = JointProb_input2/(JointProb_input1 + JointProb_input2) ; % Normalization     
          LogLikelihood(q) = real(log(JointProb_input1/JointProb_input2));
       end 
    end
    disp('Log Likelihood Ratio');
    disp(LogLikelihood);
