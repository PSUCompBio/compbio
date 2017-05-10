SET(Input_Output_src 
    Input/fe_text2matrix.cpp
    Input/fe_mainRead.cpp
    Input/Mesh.cpp
    Input/Materials.cpp
    Input/BC.cpp
    Input/Constraint.cpp
    Output/fe_vtu.cpp
    Output/fe_vector2text.cpp
    Output/fe_simulationTime.cpp
    Output/fe_matrix2text.cpp
    Output/fe_display.cpp
    Output/fe_energyWrite.cpp
     )
ADD_LIBRARY(Input-Output OBJECT ${Input_Output_src})
