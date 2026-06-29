well hello there, this is whats in this file:

this stuff is mainly based off of 
https://github.com/LiuXiangMath/Protein-Function-Prediction

and is using the version from June 29, 2026, with already the src/checkpoint things downloaded and such
# Folder Structure

nan_main/
% obv main folder thats gonna hold everything together

Protein-Function-Prediction/
% holds a clone of the LiuXiang repo (as of 06/29/26) with interproscan and the model components 

my_interproscan
% holds the stuff for interproscan

chonker_data_condensed.csv
% all proteins with annotations on 

final_output.csv
% raw combined things from the HPCC prediction runs
% around 550k rows with columns 'uniprot_id', 'sequence', and'predictions'