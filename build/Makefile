.NOTPARALLEL:
.PHONY : build build_debug run

.DEFAULT_GOAL := build

SHELL := /bin/bash

CURR_DIR := $(shell realpath ./)

CMAKE_JPEG_FLAGS := -DENABLE_SHARED=OFF -DWITH_SIMD=OFF

CXX_FLAGS := -std=c++17
DEBUG_FLAGS := -g -O0
RELEASE_FLAGS := -O3

LTO_FLAGS := -flto

MPKZEROCOST_CC := /usr/bin/clang-11
MPKZEROCOST_CXX := /usr/bin/clang++-11
MPKZEROCOST_AR := /usr/bin/llvm-ar-11
MPKZEROCOST_RANLIB := /usr/bin/llvm-ranlib-11
MPKZEROCOST_FLAGS := -fno-asm -fno-asm-blocks $(LTO_FLAGS) -fsanitize=safe-stack -fsanitize=cfi-icall -fsanitize-cfi-icall-generalize-pointers -fsanitize-cfi-canonical-jump-tables -ftrivial-auto-var-init=zero -enable-trivial-auto-var-init-zero-knowing-it-will-be-removed-from-clang

WASM_CC := /opt/wasi-sdk/bin/clang
WASM_CXX := /opt/wasi-sdk/bin/clang++
LUCET_SRC := $(shell realpath ../../lucet_sandbox_compiler)
LUCET := $(LUCET_SRC)/target/release/lucetc
LUCET_WASI_DIR := $(LUCET_SRC)/lucet-wasi
LUCET_WASI := $(LUCET_SRC)/target/release/lucet-wasi-wasmsbx

###########################`###################

stock_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS)

stock_build/libturbojpeg.a: stock_build
	cd $(shell dirname $@) && $(MAKE)

stock_build/example_image_change_quality: stock_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build/example_async_parse: stock_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build/all: stock_build/example_image_change_quality stock_build/example_async_parse

##############################################

stock_build_debug:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Debug $(CMAKE_JPEG_FLAGS)

stock_build_debug/libturbojpeg.a: stock_build_debug
	cd $(shell dirname $@) && $(MAKE)

stock_build_debug/example_image_change_quality: stock_build_debug/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build_debug/example_async_parse: stock_build_debug/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build_debug/all: stock_build_debug/example_image_change_quality stock_build_debug/example_async_parse

###########################`###################

stocklto_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_FLAGS="$(LTO_FLAGS)" \
		-DCMAKE_LD_FLAGS="$(LTO_FLAGS)"

stocklto_build/libturbojpeg.a: stocklto_build
	cd $(shell dirname $@) && $(MAKE)

stocklto_build/example_image_change_quality: stocklto_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(LTO_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stocklto_build/example_async_parse: stocklto_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(LTO_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

stock_build_lto/all: stock_build_lto/example_image_change_quality stock_build_lto/example_async_parse

##############################################

mpkzerocost_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(MPKZEROCOST_CC) \
		-DCMAKE_CXX_COMPILER=$(MPKZEROCOST_CXX) \
		-DCMAKE_AR=$(MPKZEROCOST_AR) \
		-DCMAKE_RANLIB=$(MPKZEROCOST_RANLIB) \
		-DCMAKE_C_FLAGS="$(MPKZEROCOST_FLAGS)" \
		-DCMAKE_LD_FLAGS="$(MPKZEROCOST_FLAGS)"

mpkzerocost_build/libturbojpeg.a: mpkzerocost_build
	cd $(shell dirname $@) && $(MAKE)

mpkzerocost_build/example_image_change_quality: stock_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(MPKZEROCOST_CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(MPKZEROCOST_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

mpkzerocost_build/example_async_parse: mpkzerocost_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(MPKZEROCOST_CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) $(MPKZEROCOST_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

mpkzerocost_build/all: mpkzerocost_build/example_image_change_quality mpkzerocost_build/example_async_parse

##############################################

mpkzerocost_build_debug:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Debug $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(MPKZEROCOST_CC) \
		-DCMAKE_CXX_COMPILER=$(MPKZEROCOST_CXX) \
		-DCMAKE_AR=$(MPKZEROCOST_AR) \
		-DCMAKE_RANLIB=$(MPKZEROCOST_RANLIB) \
		-DCMAKE_C_FLAGS="$(MPKZEROCOST_FLAGS)" \
		-DCMAKE_LD_FLAGS="$(MPKZEROCOST_FLAGS)"

mpkzerocost_build_debug/libturbojpeg.a: mpkzerocost_build_debug
	cd $(shell dirname $@) && $(MAKE)

mpkzerocost_build_debug/example_image_change_quality: mpkzerocost_build_debug/libturbojpeg.a ../example_image_change_quality.cpp
	$(MPKZEROCOST_CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) $(MPKZEROCOST_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@

mpkzerocost_build_debug/example_async_parse: mpkzerocost_build_debug/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(MPKZEROCOST_CXX) $(CXX_FLAGS) $(DEBUG_FLAGS) $(MPKZEROCOST_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@

mpkzerocost_build_debug/all: mpkzerocost_build_debug/example_image_change_quality mpkzerocost_build_debug/example_async_parse

##############################################

wasm_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(WASM_CC) \
		-DCMAKE_CXX_COMPILER=$(WASM_CXX) \
		-DCMAKE_C_FLAGS="-I $(shell realpath ./wasm_extra)"

wasm_build/libturbojpeg.a: wasm_build
	cd $(shell dirname $@) && $(MAKE)

wasm_build/example_image_change_quality: wasm_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(WASM_CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) -I $(shell dirname $@) ../example_image_change_quality.cpp $(shell dirname $@)/libturbojpeg.a -o $@.wasm
	$(LUCET) --bindings ${LUCET_WASI_DIR}/bindings.json \
		--guard-size "4GiB" --min-reserved-size "4GiB" --max-reserved-size "4GiB" \
		$@.wasm -o $@

wasm_build/example_async_parse: wasm_build/libturbojpeg.a $(CURR_DIR)/../example_async_parse.cpp
	$(WASM_CXX) $(CXX_FLAGS) $(RELEASE_FLAGS) \
		-I $(CURR_DIR)/../../rlbox_sandboxing_api/code/include \
		-I $(shell dirname $@) \
		$(CURR_DIR)/../example_async_parse.cpp $(shell dirname $@)/libturbojpeg.a -o $@.wasm
	$(LUCET) --bindings ${LUCET_WASI_DIR}/bindings.json \
		--guard-size "4GiB" --min-reserved-size "4GiB" --max-reserved-size "4GiB" \
		$@.wasm -o $@

wasm_build/all: wasm_build/example_image_change_quality wasm_build/example_async_parse

##############################################

build: stock_build/all stocklto_build/all wasm_build/all mpkzerocost_build/all

build_debug: stock_build_debug/all mpkzerocost_build_debug/all

run:
	@echo "Stock"                 && taskset -c 1 stock_build/example_image_change_quality
	@echo "StockLTO"              && taskset -c 1 stocklto_build/example_image_change_quality
	@echo "Wasm"                  && taskset -c 1 $(LUCET_WASI) wasm_build/example_image_change_quality
	@echo "MPKZerocost"           && taskset -c 1 mpkzerocost_build/example_image_change_quality

clean:
	-rm -rf ./*_build
	-rm -rf ./*_build_debug
