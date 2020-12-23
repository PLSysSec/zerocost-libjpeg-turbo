.NOTPARALLEL:
.PHONY : build build_debug run

.DEFAULT_GOAL := build

SHELL := /bin/bash

CURR_DIR := $(shell realpath ./)

CMAKE_JPEG_FLAGS := -DENABLE_SHARED=OFF -DWITH_SIMD=OFF
JPEG_C_FLAGS :=

CXX_FLAGS := -std=c++17
DEBUG_FLAGS := -g -O0
RELEASE_FLAGS := -O3

LTO_FLAGS := -flto
BUILD32_FLAGS := -m32

ZEROCOST_COMPILER_PATH := $(shell realpath ../../ffbuilds/zerocost_llvm_install)
# ZEROCOST_COMPILER_PATH := $(shell realpath /mnt/sata/ffbuilds/zerocost_llvm_install)
COMPILER_LIB_PATH := $(ZEROCOST_COMPILER_PATH)/lib/clang/12.0.0
ZEROCOST_CC := $(ZEROCOST_COMPILER_PATH)/bin/clang
ZEROCOST_CXX := $(ZEROCOST_COMPILER_PATH)/bin/clang++
ZEROCOST_AR := $(ZEROCOST_COMPILER_PATH)/bin/llvm-ar
ZEROCOST_RANLIB := $(ZEROCOST_COMPILER_PATH)/bin/llvm-ranlib
SEGMENT_SFI_FLAGS := -m32
ZEROCOST_COMPILE_FLAGS :=                                                            \
  -fno-asm -fno-asm-blocks                                                           \
  -fsanitize=safe-stack -fstack-clash-protection                                     \
  -flto                                                                              \
  -fsanitize=cfi-icall -fsanitize-cfi-canonical-jump-tables -fsanitize-cfi-cross-dso \
  -ftrivial-auto-var-init=zero -enable-trivial-auto-var-init-zero-knowing-it-will-be-removed-from-clang
ZEROCOST_LINK_FLAGS :=                                                               \
  -fno-asm -fno-asm-blocks                                                           \
  -fsanitize=safe-stack -fstack-clash-protection                                     \
  -flto -fuse-ld=lld                                                                 \
  -fsanitize=cfi-icall -fsanitize-cfi-canonical-jump-tables -fsanitize-cfi-cross-dso \
  ${COMPILER_LIB_PATH}/lib/linux/libclang_rt.cfi-i386.a
ZEROCOST_COMBINED_FLAGS :=                                                           \
  -fno-asm -fno-asm-blocks                                                           \
  -fsanitize=safe-stack -fstack-clash-protection                                     \
  -flto -fuse-ld=lld                                                                 \
  -fsanitize=cfi-icall -fsanitize-cfi-canonical-jump-tables -fsanitize-cfi-cross-dso \
  -ftrivial-auto-var-init=zero -enable-trivial-auto-var-init-zero-knowing-it-will-be-removed-from-clang \
  ${COMPILER_LIB_PATH}/lib/linux/libclang_rt.cfi-i386.a


WASM_CC := /opt/wasi-sdk/bin/clang
WASM_CXX := /opt/wasi-sdk/bin/clang++
LUCET_SRC := $(shell realpath ../../lucet_sandbox_compiler)
LUCET := $(LUCET_SRC)/target/release/lucetc
LUCET_WASI_DIR := $(LUCET_SRC)/lucet-wasi
LUCET_WASI := $(LUCET_SRC)/target/release/lucet-wasi-wasmsbx

###########################`###################

stock_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="$(JPEG_C_FLAGS)" $(CMAKE_JPEG_FLAGS)

stock_build/libturbojpeg.a: stock_build
	cd $(shell dirname $@) && $(MAKE) -j8

stock_build/example_image_change_quality: stock_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build/example_async_parse: stock_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build/all: stock_build/example_image_change_quality

##############################################

stock_build_debug:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS="$(JPEG_C_FLAGS)" $(CMAKE_JPEG_FLAGS)

stock_build_debug/libturbojpeg.a: stock_build_debug
	cd $(shell dirname $@) && $(MAKE) -j8

stock_build_debug/example_image_change_quality: stock_build_debug/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build_debug/example_async_parse: stock_build_debug/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build_debug/all: stock_build_debug/example_image_change_quality

###########################`###################

stocklto_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_FLAGS="$(LTO_FLAGS) $(JPEG_C_FLAGS)"

stocklto_build/libturbojpeg.a: stocklto_build
	cd $(shell dirname $@) && $(MAKE) -j8

stocklto_build/example_image_change_quality: stocklto_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(LTO_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stocklto_build/example_async_parse: stocklto_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(LTO_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stocklto_build/all: stocklto_build/example_image_change_quality

###########################`###################

stock32_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_FLAGS="$(JPEG_C_FLAGS) $(BUILD32_FLAGS)" $(CMAKE_JPEG_FLAGS)

stock32_build/libturbojpeg.a: stock32_build
	cd $(shell dirname $@) && $(MAKE) -j8

stock32_build/example_image_change_quality: stock32_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(BUILD32_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock32_build/example_async_parse: stock32_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(BUILD32_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock32_build/all: stock32_build/example_image_change_quality

##############################################

stock32_build_debug:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Debug -DCMAKE_C_FLAGS="$(JPEG_C_FLAGS) $(BUILD32_FLAGS)" $(CMAKE_JPEG_FLAGS)

stock32_build_debug/libturbojpeg.a: stock32_build_debug
	cd $(shell dirname $@) && $(MAKE) -j8

stock32_build_debug/example_image_change_quality: stock32_build_debug/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) $(BUILD32_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock32_build_debug/example_async_parse: stock32_build_debug/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) $(BUILD32_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock32_build_debug/all: stock32_build_debug/example_image_change_quality

##############################################

segmentsfizerocost_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(ZEROCOST_CC) \
		-DCMAKE_CXX_COMPILER=$(ZEROCOST_CXX) \
		-DCMAKE_AR=$(ZEROCOST_AR) \
		-DCMAKE_RANLIB=$(ZEROCOST_RANLIB) \
		-DCMAKE_C_FLAGS="$(ZEROCOST_COMBINED_FLAGS) $(JPEG_C_FLAGS) $(SEGMENT_SFI_FLAGS)"

segmentsfizerocost_build/libturbojpeg.a: segmentsfizerocost_build
	cd $(shell dirname $@) && $(MAKE) -j8

segmentsfizerocost_build/example_image_change_quality: segmentsfizerocost_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(ZEROCOST_CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(ZEROCOST_COMBINED_FLAGS) $(SEGMENT_SFI_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

segmentsfizerocost_build/example_async_parse: segmentsfizerocost_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp -o $@

segmentsfizerocost_build/all: segmentsfizerocost_build/example_image_change_quality

##############################################

segmentsfizerocost_build_debug:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Debug $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(ZEROCOST_CC) \
		-DCMAKE_CXX_COMPILER=$(ZEROCOST_CXX) \
		-DCMAKE_AR=$(ZEROCOST_AR) \
		-DCMAKE_RANLIB=$(ZEROCOST_RANLIB) \
		-DCMAKE_C_FLAGS="$(ZEROCOST_COMBINED_FLAGS) $(JPEG_C_FLAGS) $(SEGMENT_SFI_FLAGS)"

segmentsfizerocost_build_debug/libturbojpeg.a: segmentsfizerocost_build_debug
	cd $(shell dirname $@) && $(MAKE) -j8

segmentsfizerocost_build_debug/example_image_change_quality: segmentsfizerocost_build_debug/libturbojpeg.a ../example_image_change_quality.cpp
	$(ZEROCOST_CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) $(ZEROCOST_COMBINED_FLAGS) $(SEGMENT_SFI_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

segmentsfizerocost_build_debug/example_async_parse: segmentsfizerocost_build_debug/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

segmentsfizerocost_build_debug/all: segmentsfizerocost_build_debug/example_image_change_quality

##############################################

wasm_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(WASM_CC) \
		-DCMAKE_CXX_COMPILER=$(WASM_CXX) \
		-DCMAKE_C_FLAGS="-I $(shell realpath ./wasm_extra) $(JPEG_C_FLAGS)"

wasm_build/libturbojpeg.a: wasm_build
	cd $(shell dirname $@) && $(MAKE) -j8

wasm_build/example_image_change_quality: wasm_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(WASM_CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@.wasm
	$(LUCET) --bindings ${LUCET_WASI_DIR}/bindings.json \
		--guard-size "4GiB" --min-reserved-size "4GiB" --max-reserved-size "4GiB" \
		$@.wasm -o $@

wasm_build/example_async_parse: wasm_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		-o $@

wasm_build/all: wasm_build/example_image_change_quality

##############################################

build: stock_build/all stocklto_build/all stock32_build/all segmentsfizerocost_build/all wasm_build/all

build_debug: stock_build_debug/all stock32_build_debug/all segmentsfizerocost_build_debug/all

run:
	@echo "Stock"                 && taskset -c 1 stock_build/example_image_change_quality
	@echo "StockLTO"              && taskset -c 1 stocklto_build/example_image_change_quality
	@echo "Stock32LTO"            && taskset -c 1 stock32_build/example_image_change_quality
	@echo "SegmentSfiZerocost"    && taskset -c 1 segmentsfizerocost_build/example_image_change_quality
	@echo "Wasm"                  && taskset -c 1 $(LUCET_WASI) wasm_build/example_image_change_quality

clean:
	-rm -rf ./*_build
	-rm -rf ./*_build_debug
