.NOTPARALLEL:
.PHONY : build build_debug run

.DEFAULT_GOAL := build

SHELL := /bin/bash

CURR_DIR := $(shell realpath ./)

CMAKE_JPEG_FLAGS := -DENABLE_SHARED=OFF -DWITH_SIMD=OFF

LTO_FLAGS := -flto

MPKZEROCOST_CC := clang
MPKZEROCOST_CXX := clang++
MPKZEROCOST_AR := /usr/bin/llvm-ar
MPKZEROCOST_RANLIB := /usr/bin/llvm-ranlib
MPKZEROCOST_FLAGS := -fno-asm -fno-asm-blocks -flto -fsanitize=cfi-icall -fsanitize-cfi-icall-generalize-pointers -fno-sanitize-cfi-cross-dso -ftrivial-auto-var-init=zero -enable-trivial-auto-var-init-zero-knowing-it-will-be-removed-from-clang

SAFESTACK_FLAGS=-fsanitize=safe-stack

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
	$(CXX) ../example_image_change_quality.cpp -O3 -I $(shell dirname $@) -L $(shell dirname $@) -l:libturbojpeg.a -o $@

##############################################

stock_build_debug:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Debug $(CMAKE_JPEG_FLAGS)

stock_build_debug/libturbojpeg.a: stock_build_debug
	cd $(shell dirname $@) && $(MAKE)

stock_build_debug/example_image_change_quality: stock_build_debug/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) ../example_image_change_quality.cpp -g -O0 -I $(shell dirname $@) -L $(shell dirname $@) -l:libturbojpeg.a -o $@

###########################`###################

stocklto_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_FLAGS="$(LTO_FLAGS)" \
		-DCMAKE_LD_FLAGS="$(LTO_FLAGS)"


stocklto_build/libturbojpeg.a: stocklto_build
	cd $(shell dirname $@) && $(MAKE)

stocklto_build/example_image_change_quality: stocklto_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(CXX) ../example_image_change_quality.cpp $(LTO_FLAGS) -O3 -I $(shell dirname $@) -L $(shell dirname $@) -l:libturbojpeg.a -o $@

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

mpkzerocost_build/example_image_change_quality: mpkzerocost_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(MPKZEROCOST_CXX) ../example_image_change_quality.cpp $(MPKZEROCOST_FLAGS) -O3 -I $(shell dirname $@) -L $(shell dirname $@) -l:libturbojpeg.a -o $@

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
	$(MPKZEROCOST_CXX) ../example_image_change_quality.cpp $(MPKZEROCOST_FLAGS) -g -O0 -I $(shell dirname $@) -L $(shell dirname $@) -l:libturbojpeg.a -o $@

##############################################

mpkzerocostsafestack_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(MPKZEROCOST_CC) \
		-DCMAKE_CXX_COMPILER=$(MPKZEROCOST_CXX) \
		-DCMAKE_AR=$(MPKZEROCOST_AR) \
		-DCMAKE_RANLIB=$(MPKZEROCOST_RANLIB) \
		-DCMAKE_C_FLAGS="$(MPKZEROCOST_FLAGS) $(SAFESTACK_FLAGS)" \
		-DCMAKE_LD_FLAGS="$(MPKZEROCOST_FLAGS) $(SAFESTACK_FLAGS)"

mpkzerocostsafestack_build/libturbojpeg.a: mpkzerocostsafestack_build
	cd $(shell dirname $@) && $(MAKE)

mpkzerocostsafestack_build/example_image_change_quality: mpkzerocostsafestack_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(MPKZEROCOST_CXX) ../example_image_change_quality.cpp $(MPKZEROCOST_FLAGS) $(SAFESTACK_FLAGS) -O3 -I $(shell dirname $@) -L $(shell dirname $@) -l:libturbojpeg.a -o $@

##############################################

wasm_build:
	cmake -S ../ -B $@ -DCMAKE_BUILD_TYPE=Release $(CMAKE_JPEG_FLAGS) \
		-DCMAKE_C_COMPILER=$(WASM_CC) \
		-DCMAKE_CXX_COMPILER=$(WASM_CXX) \
		-DCMAKE_C_FLAGS="-I $(shell realpath ./wasm_extra)"

wasm_build/libturbojpeg.a: wasm_build
	cd $(shell dirname $@) && $(MAKE)

wasm_build/example_image_change_quality: wasm_build/libturbojpeg.a ../example_image_change_quality.cpp
	$(WASM_CXX) ../example_image_change_quality.cpp -O3 -I $(shell dirname $@) -L $(shell dirname $@) -lturbojpeg -o $@.wasm
	$(LUCET) --bindings ${LUCET_WASI_DIR}/bindings.json \
		--guard-size "4GiB" --min-reserved-size "4GiB" --max-reserved-size "4GiB" \
		$@.wasm -o $@

##############################################

build: stock_build/example_image_change_quality stocklto_build/example_image_change_quality wasm_build/example_image_change_quality mpkzerocost_build/example_image_change_quality mpkzerocostsafestack_build/example_image_change_quality

build_debug: stock_build_debug/example_image_change_quality mpkzerocost_build_debug/example_image_change_quality

run:
	@echo "Stock"                 && taskset -c 1 stock_build/example_image_change_quality
	@echo "StockLTO"              && taskset -c 1 stocklto_build/example_image_change_quality
	@echo "Wasm"                  && taskset -c 1 $(LUCET_WASI) wasm_build/example_image_change_quality
	@echo "MPKZerocost"           && taskset -c 1 mpkzerocost_build/example_image_change_quality
	@echo "MPKZerocost+Safestack" && taskset -c 1 mpkzerocostsafestack_build/example_image_change_quality

clean:
	-rm -rf ./*_build
	-rm -rf ./*_build_debug
