
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)
PROJECT(ppppppppp)
#SET(CMAKE_CXX_COMPILER "/opt/rh/devtoolset-7/root/usr/bin/g++")
SET(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR} )

FIND_PACKAGE(Nlopt)
FIND_PACKAGE(ROOT 6.14.202 REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -O3 -fopenmp  -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-int-in-bool-context -Wno-misleading-indentation")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra -Werror -O3")

SET( INCLUDE_DIRS 
    ${INCLUDE_DIRS}
    ${Nlopt_INCLUDE_DIR} 
    ${ROOT_INCLUDE_DIR} 
    ${PROJECT_BINARY_DIR}  
    ${PROJECT_SOURCE_DIR}
    $ENV{EIGEN_INCLUDE_DIR}
)

INCLUDE_DIRECTORIES( ${INCLUDE_DIRS} )

#set libraries
ADD_LIBRARY(massShape           SHARED massShape.cc           )
ADD_LIBRARY(angularDependence   SHARED angularDependence.cc   )
ADD_LIBRARY(amplitude           SHARED amplitude.cc           )
ADD_LIBRARY(integrator          SHARED integrator.cc          )
ADD_LIBRARY(generator           SHARED generator.cc           )
ADD_LIBRARY(kinematicSignature  SHARED kinematicSignature.cc  )
ADD_LIBRARY(modelAmplitude      SHARED modelAmplitude.cc      )
ADD_LIBRARY(modelGenerator      SHARED modelGenerator.cc      )
ADD_LIBRARY(logLikelihood       SHARED logLikelihood.cc       )
ADD_LIBRARY(efficiencyFunction  SHARED efficiencyFunction.cc  )
ADD_LIBRARY(BELLE_efficiency    SHARED BELLE_efficiency.cc    )
ADD_LIBRARY(getBELLEdata        SHARED getBELLEdata.cc        )
ADD_LIBRARY(phaseSpaceFunctions SHARED phaseSpaceFunctions.cc )
ADD_LIBRARY(ROOT_fit_wrapper    SHARED ROOT_fit_wrapper.cc    )

#set executables
ADD_EXECUTABLE(sandbox.exe sandbox.cxx)

# link libraries
TARGET_LINK_LIBRARIES(sandbox.exe massShape          )
TARGET_LINK_LIBRARIES(sandbox.exe angularDependence  )
TARGET_LINK_LIBRARIES(sandbox.exe amplitude          )
TARGET_LINK_LIBRARIES(sandbox.exe generator          )
TARGET_LINK_LIBRARIES(sandbox.exe integrator         )
TARGET_LINK_LIBRARIES(sandbox.exe modelAmplitude     )
TARGET_LINK_LIBRARIES(sandbox.exe modelGenerator     )
TARGET_LINK_LIBRARIES(sandbox.exe logLikelihood      )
TARGET_LINK_LIBRARIES(sandbox.exe efficiencyFunction )
TARGET_LINK_LIBRARIES(sandbox.exe getBELLEdata       )
TARGET_LINK_LIBRARIES(sandbox.exe ROOT_fit_wrapper   )

TARGET_LINK_LIBRARIES(efficiencyFunction BELLE_efficiency)

TARGET_LINK_LIBRARIES(amplitude massShape          )
TARGET_LINK_LIBRARIES(amplitude angularDependence  )
TARGET_LINK_LIBRARIES(amplitude kinematicSignature )
TARGET_LINK_LIBRARIES(amplitude kinematicSignature )

TARGET_LINK_LIBRARIES(ROOT_fit_wrapper amplitude    )
TARGET_LINK_LIBRARIES(ROOT_fit_wrapper ${ROOT_LIBS} )

TARGET_LINK_LIBRARIES(generator kinematicSignature )

TARGET_LINK_LIBRARIES(angularDependence kinematicSignature )

TARGET_LINK_LIBRARIES(modelAmplitude amplitude )

TARGET_LINK_LIBRARIES(modelGenerator modelAmplitude )
TARGET_LINK_LIBRARIES(modelGenerator generator      )

TARGET_LINK_LIBRARIES(logLikelihood amplitude     )
TARGET_LINK_LIBRARIES(logLikelihood integrator    )
TARGET_LINK_LIBRARIES(logLikelihood ${Nlopt_LIBS} )
#TARGET_LINK_LIBRARIES(logLikelihood ${ROOT_LIBS}  )

TARGET_LINK_LIBRARIES(integrator amplitude          )
TARGET_LINK_LIBRARIES(integrator generator          )
TARGET_LINK_LIBRARIES(integrator efficiencyFunction )

TARGET_LINK_LIBRARIES(getBELLEdata ${ROOT_LIBS}  )

