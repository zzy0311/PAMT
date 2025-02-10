

ifeq ($(OS),Windows_NT)
    IS_WINDOWS := yes
    RM := del *.o ptmt.exe 
else
    IS_WINDOWS := no
    RM := rm -f *.o ptmt 
endif


CXXFLAGS = -std=c++17 -fopenmp

ptmt: main.o partition.o ptmt.o 
	g++ $(CXXFLAGS) main.o partition.o ptmt.o -o PTMT

main.o: main.cpp
	g++ $(CXXFLAGS) -c main.cpp

partition.o: partition.cpp
	g++ $(CXXFLAGS) -c partition.cpp

ptmt.o: ptmt.cpp
	g++ $(CXXFLAGS) -c ptmt.cpp

clean:
	$(RM)
	