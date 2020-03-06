CXX		  := g++
#CXX_WIN	  := x86_64-w64-mingw32-g++
CXX_FLAGS := -Wall -Wextra -std=c++17 -ggdb -static

BIN		:= bin
SRC		:= src
INCLUDE	:= $FFTW_INCLUDE
LIB		:= $FFTW_LIB

LIBRARIES	:= -lfftw3
EXECUTABLE	:= main
WIN_EXECUTABLE := main_win.exe


all: $(BIN)/$(EXECUTABLE) #$(BIN)/$(WIN_EXECUTABLE)

run: clean all
	clear
	./$(BIN)/$(EXECUTABLE)

$(BIN)/$(EXECUTABLE): $(SRC)/*.cpp
	$(CXX) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

#$(BIN)/$(WIN_EXECUTABLE): $(SRC)/*.cpp
#	$(CXX_WIN) $(CXX_FLAGS) -I$(INCLUDE) -L$(LIB) $^ -o $@ $(LIBRARIES)

clean:
	-rm $(BIN)/*
