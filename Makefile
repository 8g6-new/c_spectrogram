CC       = gcc
CFLAGS   = -Wall -Wextra -O3 -march=native -mtune=native -ffast-math \
           -funroll-loops -fpeel-loops -ftracer -ftree-vectorize \
           -ftree-loop-vectorize -fopt-info-vec-optimized -fopenmp \
           -mavx -mavx2 -mfma -msse4.2 -DDEBUG 

# Debug build flags (for debugging)
CFLAGS_DEBUG = -g -O0 -DDEBUG

# Link-time optimization
CFLAGS  += -flto=auto -fuse-linker-plugin

# Profile-guided optimization flags (uncomment when using PGO)
# CFLAGS  += -fprofile-generate  # First pass
# CFLAGS  += -fprofile-use       # Second pass

# Linker flags
LDFLAGS  =  -lm -lfftw3 -lfftw3f  -lsndfile -lpng \
		   -mavx -mavx2 -mfma -g

# Target names
TARGET   = main
SOURCES  = main.c src/libheatmap/heatmap.c src/png_tools/png_tools.c src/stft/audio_tools.c
HEADERS  = src/libheatmap/heatmap.h src/png_tools/png_tools.h src/stft/audio_tools.h
OBJECTS  = $(SOURCES:.c=.o)

# Default rule for Release Build
all: $(TARGET)

# Rule to build the target executable in release mode
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# Rule to build object files in release mode
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Debug build (using debug flags)
debug: CFLAGS = $(CFLAGS_DEBUG)
debug: clean $(TARGET)
	@echo "Built in Debug Mode"
	@./$(TARGET)

run:
	./$(TARGET) ./tests/files/black_woodpecker.wav black_woodpecker 2048 128 hann 512 128

clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all run clean debug
