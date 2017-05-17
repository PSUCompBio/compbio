SET(Math_src 
    GeneralMath/fe_newtonRhapson.cpp
    GeneralMath/fe_invMatrix.cpp
    GeneralMath/fe_function.cpp
    GeneralMath/fe_find.cpp
    GeneralMath/fe_detMatrix.cpp
    GeneralMath/fe_concatenate_vector2matrix.cpp
    GeneralMath/fe_calVolume.cpp
    GeneralMath/fe_calArea_4.cpp
    GeneralMath/fe_crossProduct.cpp
    GeneralMath/fe_sort.cpp
     )
ADD_LIBRARY(Math OBJECT ${Math_src})
