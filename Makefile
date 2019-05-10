CXX	:=	g++
CXXFLAGS := -Wall -Wextra -Werror -lGL -lGLU -lglut
TARGET	:=	prog

only: main.cpp
	$(CXX) $(CXXFLAGS) main.cpp -o exec

all: cylinder.o duoCylinder.o holder.o main.o mesh.o rectangular.o specialOval.o
	$(CXX) $(CXXFLAGS) cylinder.o duoCylinder.o holder.o main.o mesh.o rectangular.o specialOval.o -o exec

cylinder: cylinder.cpp mesh.h
	$(CXX) -c cylinder.cpp -o cylinder.o

duoCylinder: duoCylinder.cpp mesh.h
	$(CXX) -c duoCylinder.cpp -o duoCylinder.o

holder: holder.cpp mesh.h
	$(CXX) -c holder.cpp -o holder.o

main: main.cpp mesh.h
	$(CXX) -c main.cpp -o main.o

mesh: mesh.cpp mesh.h
	$(CXX) -c mesh.cpp -o mesh.o

rectangular: rectangular.cpp mesh.h
	$(CXX) -c rectangular.cpp -o rectangular.o

specialOval: specialOval.cpp mesh.h
	$(CXX) -c specialOval.cpp -o specialOval.o

clean:
	rm *.o





