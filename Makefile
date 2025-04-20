# Compiler Settings
CC       = gcc
# beyond O2 its found that NAN == NAN (UB)

#_______________________General optimization flags__________________________________
OPTFLAGS+=-ffast-math#                   => Enables aggressive optimizations beyond -O2, including unsafe floating-point optimizations
OPTFLAGS+=-march=native#            => Generates code optimized for the host CPU
OPTFLAGS+=-mtune=native#            => Optimizes code scheduling for the host CPU
OPTFLAGS+=-funroll-loops#           => Unrolls loops to reduce branch overhead
OPTFLAGS+=-fpeel-loops#             => Extracts loop iterations that always execute to optimize performance
OPTFLAGS+=-ftracer#                 => Improves branch prediction and inlining
OPTFLAGS+=-flto#                    => Enables Link-Time Optimization (LTO) for cross-file optimization
OPTFLAGS+=-fuse-linker-plugin#      => Uses a linker plugin for better LTO performance
OPTFLAGS+=-MMD -MP#                 => Generates dependency files for make without including system headers
OPTFLAGS+=-floop-block#             => Optimizes loop memory access patterns for better cache performance
OPTFLAGS+=-floop-interchange#       => Reorders nested loops for better vectorization and locality
OPTFLAGS+=-floop-unroll-and-jam#    => Unrolls outer loops and fuses iterations for better performance
OPTFLAGS+=-fipa-pta#                => Enables interprocedural pointer analysis for optimization
OPTFLAGS+=-fipa-cp#                 => Performs constant propagation across functions
OPTFLAGS+=-fipa-sra#                => Optimizes function arguments and return values for efficiency
OPTFLAGS+=-fipa-icf#                => Merges identical functions to reduce code size

#_______________________Says compiler to vectorize loops_________________________________________
VECFLAGS+=-ftree-vectorize#         => Enables automatic vectorization of loops
VECFLAGS+=-ftree-loop-vectorize#    => Enables loop-based vectorization
VECFLAGS+=-fopt-info-vec-optimized# => Outputs details of vectorized loops

#_______________________SIMD Instructions that are used for vectorizing _________________________________________
VECFLAGS+=-mavx#                    => SIMD flags (4 Byte)
VECFLAGS+=-msse4.2#                 => Enables SSE4.2 instruction set (4 Byte)
VECFLAGS+=-mavx2#                   => SIMD flags (8 Byte)
VECFLAGS+=-mfma#                    => Enables Fused Multiply-Add (FMA) instructions 
VECFLAGS+=-mf16c#                   => Enables 16-bit floating-point conversion instructions
VECFLAGS+=-mabm#                    => Enables Advanced Bit Manipulation instructions
VECFLAGS+=-mf16c#                   => Enables 16-bit floating-point conversion instructions

#_______________________Debugging and safety flags__________________________________
DBGFLAGS+=-Og#                      => Optimizations suitable for debugging
DBGFLAGS+=-fno-omit-frame-pointer#  => Keeps the frame pointer for debugging
DBGFLAGS+=-fno-inline#              => Disables function inlining for better debugging
DBGFLAGS+=-fstack-protector-strong# => Adds stack protection to detect buffer overflows
DBGFLAGS+=-g#                       => Generates debugging information
DBGFLAGS+=-fsanitize=address#       => Enables AddressSanitizer for memory error detection
DBGFLAGS+=-fsanitize=leak#          => Enables leak detection
DBGFLAGS+=-fsanitize=undefined#     => Enables Undefined Behavior Sanitizer (UBSan)


LIBFLAGS = -DMINIMP3_FLOAT_OUTPUT
WARNFLAGS = -Wall -Wextra 

CFLAGS = $(WARNFLAGS) $(OPTFLAGS) $(VECFLAGS) $(LIBFLAGS)
CFLAGS_DEBUG = -g -O0 -DDEBUG
LDFLAGS = -lm -lfftw3 -lfftw3f -lsndfile -lpng -fopenmp -lopenblas 

# Directory Structure
SRCDIR     = src
SCHEMEDIR  = $(SRCDIR)/libheatmap/colorschemes
BUILDDIR   = build

# Source Files
BASE_SOURCES = main.c \
               $(wildcard $(SRCDIR)/libheatmap/*.c) \
               $(wildcard $(SRCDIR)/png_tools/*.c) \
               $(wildcard $(SRCDIR)/utils/*.c) \
               $(wildcard $(SRCDIR)/audio_tools/*.c)

# Color Scheme Sources
BUILTIN_DIR    = $(SCHEMEDIR)/builtin
OPENCV_DIR     = $(SCHEMEDIR)/opencv_like
BUILTIN_SOURCES = $(wildcard $(BUILTIN_DIR)/*.c)
OPENCV_SOURCES  = $(wildcard $(OPENCV_DIR)/*.c)

# Object Files mapped to build/
BASE_OBJECTS_BUILTIN = $(patsubst %.c,$(BUILDDIR)/builtin/%.o,$(BASE_SOURCES))
BASE_OBJECTS_OPENCV = $(patsubst %.c,$(BUILDDIR)/opencv/%.o,$(BASE_SOURCES))
BUILTIN_OBJECTS = $(patsubst %.c,$(BUILDDIR)/builtin/%.o,$(BUILTIN_SOURCES))
OPENCV_OBJECTS = $(patsubst %.c,$(BUILDDIR)/opencv/%.o,$(OPENCV_SOURCES))

# Track Last Built Target
LAST_TARGET_FILE = .last_target
ifneq ($(wildcard $(LAST_TARGET_FILE)),)
    LAST_TARGET := $(shell cat $(LAST_TARGET_FILE))
else
    LAST_TARGET := builtin
endif

# Default Target
.PHONY: all clean debug debug_opencv_like debug_builtin test opencv_like builtin run shared prep_dirs

all: $(LAST_TARGET)

# Create all needed directories
prep_dirs:
	@mkdir -p $(BUILDDIR)/builtin
	@mkdir -p $(BUILDDIR)/opencv
	@mkdir -p $(dir $(BASE_OBJECTS_BUILTIN) $(BASE_OBJECTS_OPENCV) $(BUILTIN_OBJECTS) $(OPENCV_OBJECTS))

# OpenCV Color Scheme Build - Run prep_dirs first then compile objects
opencv_like: 
	@$(MAKE) prep_dirs
	@$(MAKE) $(BASE_OBJECTS_OPENCV) $(OPENCV_OBJECTS)
	$(CC) $(CFLAGS) -DOPENCV_LIKE -o $@ $(BASE_OBJECTS_OPENCV) $(OPENCV_OBJECTS) $(LDFLAGS)
	@echo "opencv_like" > $(LAST_TARGET_FILE)
	@echo "Built with OpenCV-like color scheme"

# Builtin Color Scheme Build - Run prep_dirs first then compile objects
builtin: 
	@$(MAKE) prep_dirs
	@$(MAKE) $(BASE_OBJECTS_BUILTIN) $(BUILTIN_OBJECTS)
	$(CC) $(CFLAGS) -DBUILTIN -o $@ $(BASE_OBJECTS_BUILTIN) $(BUILTIN_OBJECTS) $(LDFLAGS)
	@echo "builtin" > $(LAST_TARGET_FILE)
	@echo "Built with Builtin color scheme"

# Debug Build with OpenCV-like Color Scheme
debug_opencv_like: CFLAGS = $(CFLAGS_DEBUG) $(LIBFLAGS)
debug_opencv_like: clean
	@$(MAKE) prep_dirs
	@$(MAKE) CFLAGS="$(CFLAGS) -DOPENCV_LIKE" $(BASE_OBJECTS_OPENCV) $(OPENCV_OBJECTS)
	$(CC) $(CFLAGS) -DOPENCV_LIKE -o opencv_like $(BASE_OBJECTS_OPENCV) $(OPENCV_OBJECTS) $(LDFLAGS)
	@echo "opencv_like" > $(LAST_TARGET_FILE)
	@echo "Debug build with OpenCV-like color scheme"

# Debug Build with Builtin Color Scheme
debug_builtin: CFLAGS = $(CFLAGS_DEBUG) $(LIBFLAGS)
debug_builtin: clean
	@$(MAKE) prep_dirs
	@$(MAKE) CFLAGS="$(CFLAGS) -DBUILTIN" $(BASE_OBJECTS_BUILTIN) $(BUILTIN_OBJECTS)
	$(CC) $(CFLAGS) -DBUILTIN -o builtin $(BASE_OBJECTS_BUILTIN) $(BUILTIN_OBJECTS) $(LDFLAGS)
	@echo "builtin" > $(LAST_TARGET_FILE)
	@echo "Debug build with Builtin color scheme"

# Shared Library Build
shared:
	@$(MAKE) prep_dirs
	@$(MAKE) $(BASE_OBJECTS_BUILTIN) $(BUILTIN_OBJECTS) $(BASE_OBJECTS_OPENCV) $(OPENCV_OBJECTS)
	$(CC) -shared -o libyourlib.so $(BASE_OBJECTS_BUILTIN) $(BUILTIN_OBJECTS) $(LDFLAGS)
	@echo "Shared library libyourlib.so built successfully."

# Compilation Rules with proper defines
$(BUILDDIR)/builtin/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -DBUILTIN -MMD -MP -c $< -o $@

$(BUILDDIR)/opencv/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -DOPENCV_LIKE -MMD -MP -c $< -o $@

# Include generated dependency files
-include $(BASE_OBJECTS_BUILTIN:.o=.d)
-include $(BASE_OBJECTS_OPENCV:.o=.d)
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
	perf stat ./$$LAST_TARGET "./tests/files/30.mp3" bird 512 128 hann 128 0 7500 64 2 3 4 "./cache/FFT"

clean:
	rm -rf $(BUILDDIR)
	rm -f builtin opencv_like main $(LAST_TARGET_FILE) libyourlib.so

