# === Toolchain ===
CXX       = g++
CXXSTD    = -std=c++20
CXXWARN   = -Wall -Wextra -Wpedantic
OPTFLAGS  = -O2
PIC       = -fPIC

# === ONNX Runtime (可覆盖) ===
ORT_PREFIX ?= /opt/onnxruntime
ORT_INC     = $(ORT_PREFIX)/include
ORT_LIBDIR  = $(ORT_PREFIX)/lib

# JSON 单头文件目录（可不设，如果系统自带）
JSON_INC ?=

# 是否启用 ONNX（源码需配合宏，见下）1
WITH_ONNX ?= 1

# === Targets ===
SRCS      = PtCalculator.cc
OBJS      = $(SRCS:.cc=.o)
TARGET_SO = libPtCalc.so

# === Includes ===
INCFLAGS  = -I.
ifeq ($(WITH_ONNX),1)
  	INCFLAGS += -I$(ORT_INC)
endif
ifneq ($(strip $(JSON_INC)),)
  	INCFLAGS += -I$(JSON_INC)
endif

# === Flags ===
CXXFLAGS  = $(CXXSTD) $(CXXWARN) $(OPTFLAGS) $(PIC) $(INCFLAGS)
LDFLAGS   = -shared
LDLIBS    =

ifeq ($(WITH_ONNX),1)
	CXXFLAGS += -DPTCALC_WITH_ONNX=1
	LDFLAGS  += -Wl,-rpath,$(ORT_LIBDIR) -L$(ORT_LIBDIR)
	LDLIBS   += -lonnxruntime
else
  	CXXFLAGS += -DPTCALC_WITH_ONNX=0
endif

# === Rules ===
all: $(TARGET_SO)

$(TARGET_SO): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)
	@echo "Built $@"

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET_SO)

print-%:
	@echo '$*=$($*)'
