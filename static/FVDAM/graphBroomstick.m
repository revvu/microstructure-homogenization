% Reevu Adakroy
% 2/8/2022
SIGMA22_collection = readmatrix("SIGMA22_collection.csv");
SIGMA22_square_hex = readmatrix("SIGMA22_square_hex.csv");
%plot(SIGMA22_collection(:,1),SIGMA22_collection(:,2:end))
plot(SIGMA22_collection(:,1),SIGMA22_collection(:,2:end),'Color',"#C0C0C0")
hold on
plot(SIGMA22_square_hex(:,1),SIGMA22_square_hex(:,2:end))
hold off
xlabel("Strain \epsilon ")
ylabel("Stress \sigma (N/m^2)")
