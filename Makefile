CXX			= g++

CXXFLAGS	= -c -g -O3 -std=c++11

INCLUDE		= -I/users/student/mr111/dhlin22/ALS/Final/

LIBRARY		= -L/users/student/mr111/dhlin22/ALS/Final/GLPK/lib -lglpk -lm -L/users/student/mr111/dhlin22/ALS/Final/LEDA/lib -lG -lL -lm

SRC			= 111062684_ALS_Final.cpp

OBJ			= ${SRC:%.cpp=%.o}

RM			= rm

EXE			= ./Power_Optimization_State_Assignment

OUT			= ./results/* ./log/*

all :: $(EXE)
object: $(SRC)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c -o $(OBJ) $(SRC)
$(EXE): object
	$(CXX) -o $(EXE) $(OBJ) $(LIBRARY)
clean:
	$(RM) -rf $(OBJ) $(EXE) $(OUT)
