# Synthetic-Vessels
Summary

Computer generated vessel trees offer an innovative approach to improving the diagnosis of cardiovascular diseases through the extension of vessel trees imaged from CT angiography. The synthetically generated vessels extend the coronary arteries from the CT images down to the arteriole level allowing for improved blood flow simulations.

In this report a method to simulate realistic coronary vascular growth within a defined 2-dimensional area has been proposed and assessed. The method aims to minimise the total volume of the network whilst abiding by certain physical principles. These include angle rules to dictate the direction of vessel growth, as well as intersection rules to prevent vessels from growing within each other. This method has been used to generate vessel trees of up to 4000 terminal segments in an average time of 2 hours. The proposed method has been assessed against other leading research within the field as well as morphometric data taken from porcine hearts. In comparison to the leading research the model fits the shape of diameter plots against bifurcation level and Strahler order. However, discrepancies within the initial vessels are apparent due to inaccuracy within the seeding data. As well as this consistently the model’s vessels are smaller than that of stated literature likely due to the assumptions made within the pressure drop across vessels within the method. 

Due to these discrepancies improvements have been suggested for future work. With these improvements to the method it will produce an accurate representation of vascular growth of the coronary arteries within a computation time feasible for clinical diagnoses.

TreeGenerationScript.m is the main scrip that the simlation is run from
