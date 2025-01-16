CC       = gcc
CFLAGS   = -Wall -Wextra -O3 -march=native -mtune=native -ffast-math \
           -funroll-loops -fpeel-loops -ftracer -ftree-vectorize \
           -ftree-loop-vectorize -fopt-info-vec-optimized -fopenmp \
           -mavx -mavx2 -mfma -msse4.2 -DDEBUG 


# Link-time optimization
CFLAGS  += -flto=auto -fuse-linker-plugin

# Profile-guided optimization flags (uncomment when using PGO)
# CFLAGS  += -fprofile-generate  # First pass
# CFLAGS  += -fprofile-use       # Second pass



# Linker flags
LDFLAGS  = -Wl,-O3 -Wl,--as-needed -Wl,--strip-all \
           -lm -lfftw3 -lfftw3f  -lsndfile -lpng \
		   -mavx -mavx2 -mfma -g

# Target names
TARGET   = main
SOURCES  = main.c src/libheatmap/heatmap.c src/png_tools/png_tools.c src/stft/audio_tools.c
HEADERS  = src/libheatmap/heatmap.h src/png_tools/png_tools.h src/stft/audio_tools.h
OBJECTS  = $(SOURCES:.c=.o)

# Default rule
all: $(TARGET)

# Rule to build the target executable
$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS) $(LDFLAGS)

# Rule to build object files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@


run:
	./$(TARGET) ./tests/files/black_woodpecker.wav c_mel_black_woodpecker.png 4096 64 hann 512
#               input_file              output_file          window hop type mel_banks

clean:
	rm -f $(TARGET) $(OBJECTS)

.PHONY: all run clean
