# Configuration of the executable
TARGET = multiPatchLD

INSTALL_PATH = $(PWD)

# Compiler configuration
CXX      = g++
CXXFLAGS = -g -pg -Wcomment 
#CXXFLAGS = -Wall -Werror -Wextra -Wshadow -O3
COMP = -c 
INC =  -I/usr/local/opt/gsl/include -I/usr/include/c++/4.2.1/
INC_LIB =   -L/usr/local/opt/gsl/lib  -lgsl -lgslcblas -lm


OBJ =   association.o \
	buildDomainsOnLDParticles.o \
	calcEnergy.o \
        calcForce.o \
	createDomains.o \
	dissociation.o \
	drawAngles.o \
	drawOrientation.o \
	eulerAngles.o \
	findRoot.o \
        forwardFluxSampling.o \
        generateAndLoadParticles.o \
        gradientOfThePotential.o \
	greensFunction.o \
        integrator.o \
	linkedCell.o \
        loadSpeciesInfo.o \
        main.o \
        measureInterParticleDistance.o \
        measureInterPatchDistance.o \
        myVector.o \
	nearestNeighborSearches.o \
        particle.o \
        quaternion.o \
	removeTorqueComponent.o \
	rotationalGreensFunction.o \
        simulationParameters.o \
        species.o \
        updatePosition.o \
	wignerSmalld.o \
        writeOutput.o 
    
default: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJ) $(INC) $(INC_LIB) -pg

# clean
clean:
	@echo "clean ..."
	@$(RM) sph *.o *.md *.vtk

