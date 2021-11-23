SP = [10 10 1 0.1 0];
MB = [50 1 1 1];
d = [5.5 1 0.05];
Inertia = sim('SC_I.slx');
plot(Inertia.I_SC, 'o')
