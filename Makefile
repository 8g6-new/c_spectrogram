# Compiler Settings
CC       = gcc
CFLAGS   = -Wall -Wextra -O3 -march=native -mtune=native -ffast-math \
           -funroll-loops -fpeel-loops -ftracer -ftree-vectorize \
           -ftree-loop-vectorize -fopt-info-vec-optimized -fopenmp \
           -mavx -mavx2 -mfma -msse4.2 -DDEBUG \
           -flto=auto -fuse-linker-plugin -DMINIMP3_FLOAT_OUTPUT -fPIC -lopenblas 
		   
CFLAGS_DEBUG = -g -O0 -DDEBUG
LDFLAGS  = -lm -lfftw3 -lfftw3f -lsndfile -lpng -g

# Directory Structure
SRCDIR    = src
SCHEMEDIR = $(SRCDIR)/libheatmap/colorschemes

# Source Files
BASE_SOURCES = main.c \
               $(SRCDIR)/libheatmap/heatmap.c \
               $(SRCDIR)/libheatmap/heatmap_tools.c \
               $(SRCDIR)/png_tools/png_tools.c \
               $(SRCDIR)/utils/ftype_detect.c \
			   $(SRCDIR)/utils/bench.c \
               $(SRCDIR)/audio_tools/audio_io.c \
               $(SRCDIR)/audio_tools/audio_visualizer.c \
               $(SRCDIR)/audio_tools/spectral_features.c

# Color Scheme Sources
BUILTIN_DIR  = $(SCHEMEDIR)/builtin
OPENCV_DIR   = $(SCHEMEDIR)/opencv_like
BUILTIN_SOURCES = $(wildcard $(BUILTIN_DIR)/*.c)
OPENCV_SOURCES  = $(wildcard $(OPENCV_DIR)/*.c)

# Object Files
BASE_OBJECTS    = $(BASE_SOURCES:.c=.o)
BUILTIN_OBJECTS = $(BUILTIN_SOURCES:.c=.o)
OPENCV_OBJECTS  = $(OPENCV_SOURCES:.c=.o)

# Track Last Built Target
LAST_TARGET_FILE = .last_target
ifneq ($(wildcard $(LAST_TARGET_FILE)),)
    LAST_TARGET := $(shell cat $(LAST_TARGET_FILE))
else
    LAST_TARGET := builtin
endif

# Default Target
.PHONY: all clean debug test opencv_like builtin run shared

all: $(LAST_TARGET)

# OpenCV Color Scheme Build
opencv_like: CFLAGS += -DOPENCV_LIKE
opencv_like: $(BASE_OBJECTS) $(OPENCV_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "opencv_like" > $(LAST_TARGET_FILE)
	@echo "Built with OpenCV-like color scheme"

# Builtin Color Scheme Build
builtin: CFLAGS += -DBUILTIN
builtin: $(BASE_OBJECTS) $(BUILTIN_OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "builtin" > $(LAST_TARGET_FILE)
	@echo "Built with Builtin color scheme"

# Shared Library Build
shared: $(BASE_OBJECTS) $(BUILTIN_OBJECTS) $(OPENCV_OBJECTS)
	$(CC) -shared -o libyourlib.so $^ $(LDFLAGS)
	@echo "Shared library libyourlib.so built successfully."

# Compilation Rule with Dependency Tracking
%.o: %.c
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

# Include generated dependency files
-include $(BASE_OBJECTS:.o=.d)
-include $(BUILTIN_OBJECTS:.o=.d)
-include $(OPENCV_OBJECTS:.o=.d)

# Debug Build
debug: CFLAGS = $(CFLAGS_DEBUG)
debug: clean $(LAST_TARGET)
	@echo "Built in Debug Mode"


run:
	@if [ ! -f "$(LAST_TARGET_FILE)" ]; then \
	  echo "No previous build found. Run 'make' first."; exit 1; \
	fi; \
	LAST_TARGET=$$(cat $(LAST_TARGET_FILE)); \
	if [ "$$LAST_TARGET" = "shared" ]; then \
	  echo "Last build was a shared library. Nothing to run."; exit 1; \
	fi; \
	if [ ! -x "$$LAST_TARGET" ]; then \
	  echo "Executable '$$LAST_TARGET' not found. Run 'make' first."; exit 1; \
	fi; \
	echo "Running $$LAST_TARGET..."; \
	./$$LAST_TARGET "./tests/files/173.mp3"  black_woodpecke 2048 128 hann 512 128 7500 64 2 2 2 "./cache/FFT"

    # <ip_filename> <op_filename> <window_size> <hop_size> <window_type> <number_of_mel_banks> <min_mel> <max_mel> <num_coff


clean:
	find $(SRCDIR) -name "*.o" -type f -delete
	find $(SRCDIR) -name "*.d" -type f -delete
	rm -f builtin opencv_like main *.o *.d $(LAST_TARGET_FILE) libyourlib.so