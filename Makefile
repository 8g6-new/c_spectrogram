# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -O3 -march=native -funroll-loops

# STFT test executable
TARGET  = main
SOURCES = main.c src/stft/audio_tools.c src/png_tools/png_tools.c src/libheatmap/heatmap.c
HEADERS = src/stft/audio_tools.h /src/png_tools/png_tools.h src/libheatmap/heatmap.h
OBJECTS = $(SOURCES:.c=.o)
LDFLAGS = -lm -lfftw3 -lfftw3f -lsndfile -lpng 

# Default rule
all: $(TARGET)

# Build the STFT test executable
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# Compile source files into object files
%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

# Run the STFT test executable
run:
	./$(TARGET) ./tests/files/wind.wav 512 128 hann 256 wind.png
#               file                  window hop type mel_banks

clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all test clean
