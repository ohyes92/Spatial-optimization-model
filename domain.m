function z=domain(x)
    circuit_second=dlmread('C:\Users\mrsol\OneDrive\Desktop\connect_150/4class_connect1501.txt');
    lst_second=dlmread('C:\Users\mrsol\OneDrive\Desktop\lst_150/exact_lst1501.txt');
    unusedland_second=dlmread('C:\Users\mrsol\OneDrive\Desktop\empty lot_150/unusedland_1501.txt');
    weight=dlmread('C:\Users\mrsol\OneDrive\Desktop\connect_150/weight11.txt');
    
    z1 = -sum(circuit_second(pop(i).Position==4).*weight(pop(i).Position==4));
    
    z2 = -sum(lst_second(pop(i).Position==4));
    
    z=[z1 z2]';

end

