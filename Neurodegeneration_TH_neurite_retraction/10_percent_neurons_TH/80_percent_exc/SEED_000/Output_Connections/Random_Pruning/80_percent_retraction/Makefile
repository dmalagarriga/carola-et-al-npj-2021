####### Compiler, tools and options
TARGET = nmm.exe
SOURCES_DIR = Sources
HEADERS_DIR = Headers
OBJECTS_DIR = Objects

#C++ = mpic++
C++ = g++

MAIN = $(SOURCES_DIR)/$(TARGET:.exe=.cpp)

# Makes a list of the object (.o) files.
OBJECT_FILES = $(OBJECTS_DIR)/Parameters.o \
               $(OBJECTS_DIR)/neuron.o \
               $(OBJECTS_DIR)/Heun.o \
               $(OBJECTS_DIR)/rhs.o 

all: $(OBJECT_FILES)
#	@$(C++) -Wall $(MAIN) $(OBJECT_FILES) -o $(TARGET)  -lm
	@$(C++) -Wall $(MAIN) $(OBJECT_FILES) -o $(TARGET) -lgsl -lgslcblas -lm 

#	./$(TARGET)

# Clean
clean:
	rm -f *.exe *~ */*~ Objects/*.o

$(OBJECTS_DIR)/Parameters.o : $(SOURCES_DIR)/Parameters.cpp
	$(C++) -c $< -o $@

$(OBJECTS_DIR)/neuron.o : $(SOURCES_DIR)/neuron.cpp
	$(C++) -c $< -o $@

$(OBJECTS_DIR)/Heun.o : $(SOURCES_DIR)/Heun.cpp
	$(C++) -c $< -o $@

	
$(OBJECTS_DIR)/rhs.o : $(SOURCES_DIR)/rhs.cpp
	$(C++) -c $< -o $@


