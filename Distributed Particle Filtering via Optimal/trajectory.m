%% generate target trajectory
function state = trajectory(STI,t,sigma_u)
%

rng(0);

x = [0, 0, 4, 13, -1, -3]'; % origin state vector

state = [x,zeros(6,STI)];
for sti = 1:STI
    x = transition(x,t,sigma_u);
    state(:,sti+1) = x;
end

rng shuffle;
