ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS      = $(shell $(ROOTSYS)/bin/root-config --libs)
ROOTGLIBS     = $(shell $(ROOTSYS)/bin/root-config --glibs)

CXX           = g++
CXXFLAGS      = -g -Wall -fPIC -Wno-deprecated

NGLIBS         = $(ROOTGLIBS) 
NGLIBS        += -lMinuit -lFoam -lASImage

CXXFLAGS      += $(ROOTCFLAGS)
CXX           += -I./
LIBS           = $(ROOTLIBS) 

GLIBS          = $(filter-out -lNew, $(NGLIBS))

CXX	      += -I./lib/
OUTLIB	      = ./lib/
OUTBIN	      = ./bin/
.SUFFIXES: .cc,.C
.PREFIXES: ./lib/

EXECUTABLE = GuariguanchiApp

#----------------------------------------------------#

main: main.cc \
	Guariguanchi.o  Setup.o  GeometryHandler.o  AnalysisHandler.o  \
	GGlobalTools.o  GUtiliratyFunctions.o  \
	GSurfaceObject.o GSurfacePlane.o  GSurfaceCylinder.o  GSurfaceCylinderSection.o  GSurfaceDisk.o  GSurfaceDiskSection.o  GSurfaceCone.o  GSurfaceConeSection.o  GSurfacePetal.o  \
	GGeoObject.o  GGeoPlane.o  GGeoCylinder.o  GGeoCylinderSection.o  GGeoDisk.o  GGeoDiskSection.o  GGeoCone.o  GGeoConeSection.o  GGeoPetal.o \
	GGeometry.o  \
	GBField.o  GBFieldConstant.o  GBFieldMultipleSteps.o  \
	GResolutionModel.o  GResolutionModelTPC.o  \
	GEfficiencyModel.o  \
	GTrackFinderAlgo.o  GTrackFinderAlgoFPCCD.o  \
	GTrajectory.o  GTrajectoryStraight.o  GTrajectoryHelix.o  GTrajectoryMultipleSteps.o  \
	GTracker.o
	$(CXX) $(CXXFLAGS)  -o $(OUTBIN)$(EXECUTABLE)  $(OUTLIB)/*.o $(GLIBS)  $<

# global tools
GGlobalTools.o:  src/GGlobalTools.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGlobalTools.o $<

GUtiliratyFunctions.o:  src/GUtiliratyFunctions.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GUtiliratyFunctions.o $<


# Surface object classes
GSurfaceObject.o: src/GSurfaceObject.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceObject.o $<

GSurfacePlane.o: src/GSurfacePlane.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceOPlane.o $<

GSurfaceCylinder.o: src/GSurfaceCylinder.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceCylinder.o $<
	
GSurfaceCylinderSection.o: src/GSurfaceCylinderSection.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceCylinderSection.o $<
	
GSurfaceDisk.o: src/GSurfaceDisk.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceDisk.o $<
	
GSurfaceDiskSection.o: src/GSurfaceDiskSection.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceDiskSection.o $<
	
GSurfaceCone.o: src/GSurfaceCone.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceCone.o $<
	
GSurfaceConeSection.o: src/GSurfaceConeSection.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfaceConeSection.o $<
	
GSurfacePetal.o: src/GSurfacePetal.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GSurfacePetal.o $<

# Geometry object classes
GGeoObject.o: src/GGeoObject.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoObject.o $<

GGeoPlane.o: src/GGeoPlane.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoPlane.o $<
	
GGeoCylinder.o: src/GGeoCylinder.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoCylinder.o $<
	
GGeoCylinderSection.o: src/GGeoCylinderSection.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoCylinderSection.o $<
	
GGeoDisk.o: src/GGeoDisk.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoDisk.o $<
	
GGeoDiskSection.o: src/GGeoDiskSection.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoDiskSection.o $<
	
GGeoCone.o: src/GGeoCone.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoCone.o $<
	
GGeoConeSection.o: src/GGeoConeSection.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoConeSection.o $<
	
GGeoPetal.o: src/GGeoPetal.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeoPetal.o $<

# Geometry classe
GGeometry.o: src/GGeometry.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GGeometry.o $<

#BField classes
GBField.o: src/GBField.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GBField.o $<
	
GBFieldConstant.o: src/GBFieldConstant.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GBFieldConstant.o $<
	
GBFieldMultipleSteps.o: src/GBFieldMultipleSteps.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GBFieldMultipleSteps.o $<

	
#Resolution model classes
GResolutionModel.o: src/GResolutionModel.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GResolutionModel.o $<

GResolutionModelTPC.o: src/GResolutionModelTPC.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GResolutionModelTPC.o $<

#Efficiency model classes
GEfficiencyModel.o: src/GEfficiencyModel.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GEfficiencyModel.o $<

#Trajectory classes
GTrajectory.o: src/GTrajectory.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GTrajectory.o $<
	
GTrajectoryStraight.o: src/GTrajectoryStraight.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GTrajectoryStraight.o $<
	
GTrajectoryHelix.o: src/GTrajectoryHelix.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GTrajectoryHelix.o $<
	
GTrajectoryMultipleSteps.o: src/GTrajectoryMultipleSteps.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GTrajectoryMultipleSteps.o $<

#Track finder classes
GTrackFinderAlgo.o: src/GTrackFinderAlgo.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GTrackFinderAlgo.o $<
	
GTrackFinderAlgoFPCCD.o: src/GTrackFinderAlgoFPCCD.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GTrackFinderAlgoFPCCD.o $<
	
# Tracker class code
GTracker.o: src/GTracker.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GTracker.o $<

	
# main guarigunchi code
Guariguanchi.o:  src/Guariguanchi.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)Guariguanchi.o $<

Setup.o:  src/Setup.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)Setup.o $<

GeometryHandler.o:  src/GeometryHandler.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)GeometryHandler.o $<

AnalysisHandler.o:  src/AnalysisHandler.cc
	$(CXX) $(CXXFLAGS) -c -I. -o $(OUTLIB)AnalysisHandler.o $<

clean:
	rm -f src/*~
	rm -f *~
	rm -f $(OUTLIB)*.o
	rm -f $(OUTBIN)$(EXECUTABLE)

#----------------------------------------------------#

